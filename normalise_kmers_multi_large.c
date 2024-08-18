#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <pthread.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <time.h>
#include <locale.h>
#include <getopt.h>
#include <stdbool.h>

//////// HELPERS

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

/// registry of allocs (malloc realloc calloc)
// output_filename
// hash_table_t
// hash_table_t->entries
// forward_ptr
// reverse_ptr
// cfg.forward_files
// cfg.reverse_files

// TODO
// make every thread independent, remove global table and sync, depth divided by cpus
//
//
//

// #define INITIAL_CAPACITY 500000003 // 7.5gb per cpu equiv genome/transcriptome size; should be a prime number; remember we do not convert to canonical kmers (as this was create for rnaseq)
#define INITIAL_CAPACITY 150000001 // 2.3gb good starting value for a transcriptome, esp if stranded
#define MAX_SEQ_LENGTH 400
#define CHUNK_SIZE 1000000
#define MAX_THREADS 256
// #define SYNC_INTERVAL 100000 // how many seqs/often to sync the kmer tables of each thread, this inits the thread-specific sync interval which increases with more data.
#define SYNC_INTERVAL 1000000
#define TABLE_LOAD_FACTOR 0.8
#define MAX_K 32
#define REPORTING_INTERVAL 60 // seconds

// lookup table, including lower case (Just In Case(TM))

static const uint8_t base_map[256] = {
    ['A'] = 0x00, ['C'] = 0x01, ['G'] = 0x02, ['T'] = 0x03
    // , ['a'] = 0x00, ['c'] = 0x01, ['g'] = 0x02, ['t'] = 0x03 // assume uppercase
};
static const char rev_base_map[4] = {'A', 'C', 'G', 'T'};

// 16 (8+4+4padding) bytes per entry
typedef struct
{
    uint64_t hash;
    int count;
} kmer_t;

typedef struct hash_table_t
{
    kmer_t *entries;
    size_t used;
    size_t capacity;
} hash_table_t;

// if we decide on partitions.
// #define NUM_PARTITIONS 4
// hash_table_t global_hash_tables[NUM_PARTITIONS];

typedef struct
{
    char *forward_ptr;
    char *reverse_ptr;
    size_t forward_size;
    size_t reverse_size;
    hash_table_t *hash_table;
    // hash_table_t *hash_tables[NUM_PARTITIONS];
    int thread_id;
    size_t processed_count;
    size_t printed_count;
    size_t skipped_count;
    size_t total_kmers;
    time_t last_report_time;
    int last_report_count;
    FILE *thread_output_forward;
    FILE *thread_output_reverse;
    char *thread_output_forward_filename;
    char *thread_output_reverse_filename;
    int thread_sync_interval;
} thread_data_t;

thread_data_t thread_data[MAX_THREADS];

struct reporting_t
{
    size_t total_processed;
    size_t total_printed;
    size_t total_skipped;
    size_t total_kmers;
    int files_processed;
} reporting;

struct config_t
{
    char **forward_files;
    int forward_file_count;
    char **reverse_files;
    int reverse_file_count;
    int ksize;
    int depth;
    float coverage;
    int verbose;
    char *filetype;
    int cpus;
    int help;
    int debug;
    int canonical;
    int memory;
    size_t initial_hash_size;
} cfg;

typedef struct
{
    char *data;
    size_t size;
} mmap_file_t;

// Function prototypes
char *read_line(char *ptr, char *buffer, int max_length);
float capacity2memory(size_t capacity);
float memoryGB2capacity(size_t memory);
mmap_file_t mmap_file(const char *filename);
void munmap_file(mmap_file_t *mf);
char *find_next_record_start(char *ptr, char *end_ptr, char stop_char);
static void replacestr(const char *line, const char *search, const char *replace);

////

void print_usage(char *program_name);
int parse_arguments(int argc, char *argv[]);
void process_forward_files(char *first_file, char **additional_files);
void process_reverse_files(char *first_file, char **additional_files);
char *create_output_filename(const char *basename, int k, int norm_depth, int thread);

////

bool is_hash_table_empty(const hash_table_t *ht);
void init_hash_table(hash_table_t *ht);
hash_table_t *copy_hash_table(const hash_table_t *source);
int store_kmer(hash_table_t *hash_table, uint64_t hash, int thread_id);
size_t expand_local_hash_table(hash_table_t *ht, size_t new_capacity, int thread_id);
void synchronise_hash_tables(hash_table_t *local_ht, int thread_id);

////

uint64_t murmurhash3(const char *key, int len, uint32_t seed);
static inline uint64_t encode_kmer_murmur(const char *seq, int k);
static inline uint64_t encode_kmer_plain(const char *seq, int k);
char decode_kmer_plain(uint64_t encoded, int k, char *kmer);
static inline uint64_t mix_bits(uint64_t hash);
static inline uint64_t kmer_hash_fnv(const char *seq, int k);

////

bool valid_dna(const char *sequence);
void reverse_complement(const char *seq, char *rev_comp, int k);
const char *get_canonical_kmer(const char *kmer, int k);

////

void process_sequence(const char *seq, hash_table_t *hash_table, int K, int NORM_DEPTH, int *seq_high_count_kmers, int *total_seq_kmers, int thread_id);
void *process_thread_chunk(void *arg);
int multithreaded_process_files(thread_data_t *thread_data, mmap_file_t *forward_mmap, mmap_file_t *reverse_mmap);

////////////////////////////////////////////////////////////////
////////// HELPER FUNCTIONS ////////////////////////////////////
////////////////////////////////////////////////////////////////

char *read_line(char *ptr, char *buffer, int max_length)
{
    int i = 0;
    while (*ptr != '\n' && *ptr != '\0' && i < max_length - 1)
    {
        buffer[i++] = *ptr++;
    }
    buffer[i] = '\0';

    if (*ptr == '\n')
    {
        ptr++; // Move past the newline character
    }

    return (*ptr == '\0') ? NULL : ptr;
}

float capacity2memory(size_t capacity)
{
    float number = (float)capacity;
    return number * 16 / 1073741824;
}

float memoryGB2capacity(size_t memory)
{
    float number = (float)memory * 1073741824;
    float per_cpu = number / 16;
    return per_cpu / cfg.cpus;
}

mmap_file_t mmap_file(const char *filename)
{
    mmap_file_t mf = {NULL, 0};
    int fd = open(filename, O_RDONLY);
    if (fd == -1)
    {
        perror("Error opening file for mmap");
        return mf;
    }

    struct stat sb;
    if (fstat(fd, &sb) == -1)
    {
        perror("Error getting file size");
        close(fd);
        return mf;
    }

    mf.data = mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
    if (mf.data == MAP_FAILED)
    {
        perror("Error mmapping the file");
        close(fd);
        mf.data = NULL;
        return mf;
    }

    mf.size = sb.st_size;
    close(fd);
    return mf;
}

void munmap_file(mmap_file_t *mf)
{
    if (mf->data != NULL)
    {
        munmap(mf->data, mf->size);
        mf->data = NULL;
        mf->size = 0;
    }
}
char *find_next_record_start(char *ptr, char *end_ptr, char stop_char)
{
    while (ptr < end_ptr && *ptr != stop_char)
    {
        ptr = strchr(ptr, '\n');
        if (ptr)
            ptr++;
        else
            break;
    }
    return (ptr < end_ptr) ? ptr : NULL;
}

static void replacestr(const char *line, const char *search, const char *replace)
{
    char *sp;
    while ((sp = strstr(line, search)) != NULL)
    {
        int search_len = strlen(search);
        int replace_len = strlen(replace);
        int tail_len = strlen(sp + search_len);
        memmove(sp + replace_len, sp + search_len, tail_len + 1);
        memcpy(sp, replace, replace_len);
    }
}

////////////////////////////////////////////////////////////////
///////////////// ARGUMENTS AND FILES //////////////////////////
////////////////////////////////////////////////////////////////

void print_usage(char *program_name)
{
    fprintf(stderr, "Usage: %s\t--forward file1 [file2+] --reverse file1 [file2+] [--ksize (int; def. 25)] [--depth|-d (int; def. 100)]\n\
    \t\t\t[--coverage|g (float 0-1; def. 0.9)] [--verbose] [--filetype|-t (fq|fa; def. fq)] [--cpu|-p (int; def 1)] [--debug|-b]\n\
    [--canonical|c] [--memory|m (int; def. 150000001)] \n",
            program_name);
}

int parse_arguments(int argc, char *argv[])
{
    cfg.coverage = 0.9;
    cfg.verbose = 0;
    cfg.debug = 0;
    cfg.filetype = "fq";
    cfg.cpus = 1;
    cfg.help = 0;
    cfg.forward_file_count = 0;
    cfg.reverse_file_count = 0;
    cfg.forward_files = NULL;
    cfg.reverse_files = NULL;
    cfg.ksize = 25;
    cfg.depth = 100;
    cfg.memory = 0;
    cfg.canonical = 0;

    static struct option long_options[] = {
        {"forward", required_argument, 0, 'f'},
        {"reverse", required_argument, 0, 'r'},
        {"ksize", required_argument, 0, 'k'},
        {"depth", required_argument, 0, 'd'},
        {"coverage", required_argument, 0, 'g'},
        {"filetype", required_argument, 0, 't'},
        {"cpu", required_argument, 0, 'p'},
        {"memory", required_argument, 0, 'm'},
        {"debug", required_argument, 0, 'b'},
        {"verbose", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        {"canonical", no_argument, 0, 'c'},
        {0, 0, 0, 0}};

    int opt;
    int option_index = 0;

    while ((opt = getopt_long(argc, argv, "f:r:k:d:g:t:p:m:b:vhc", long_options, &option_index)) != -1)
    {
        switch (opt)
        {
        case 'c':
            cfg.canonical = 1;
            break;
        case 'm':
            cfg.memory = atoi(optarg);
            if (cfg.memory < 1)
            {
                printf("Memory cannot be less than 1 Gb %'d\n", cfg.memory);
                return 0;
            }
            break;
        case 'b':
            cfg.debug = atoi(optarg);
            break;
        case 'h':
            cfg.help = 1;
            return 0;
        case 'p':
            cfg.cpus = atoi(optarg);
            break;
        case 'f':
            process_forward_files(optarg, &argv[optind]);
            break;
        case 'r':
            process_reverse_files(optarg, &argv[optind]);
            break;
        case 'k':
            cfg.ksize = atoi(optarg);
            break;
        case 'd':
            cfg.depth = atoi(optarg);
            break;
        case 'g':
            cfg.coverage = atof(optarg);
            break;
        case 'v':
            cfg.verbose = 1;
            break;
        case 't':
            cfg.filetype = optarg;
            break;
        default:
            fprintf(stderr, "Unexpected error in option processing\n");
            return 0;
        }
    }

    if (cfg.verbose)
    {
        printf("CMD: ");
        for (int i = 0; i < argc; i++)
        {
            printf("%s ", argv[i]);
        }
        printf("\n\n");
    }
    if (cfg.debug > 1)
    {
        printf("Provided forward files: ");
        for (int i = 0; i < cfg.forward_file_count; i++)
        {
            printf("%s ", cfg.forward_files[i]);
        }
        printf("\nProvided reverse files: ");
        for (int i = 0; i < cfg.reverse_file_count; i++)
        {
            printf("%s ", cfg.reverse_files[i]);
        }
        printf("\n");
    }

    if (cfg.forward_file_count == 0 || cfg.reverse_file_count == 0 || cfg.ksize <= 0 || cfg.depth <= 0 || cfg.coverage < 0 || cfg.coverage > 1)
    {
        printf("ERROR: Provided options: fwd %d rev %d ksize %d depth %d coverage %.2f\n", cfg.forward_file_count, cfg.reverse_file_count, cfg.ksize, cfg.depth, cfg.coverage);
        return 0;
    }

    if (cfg.forward_file_count != cfg.reverse_file_count)
    {
        fprintf(stderr, "Error: Number of forward and reverse files must match\n");
        return 0;
    }
    if (cfg.cpus <= 0)
    {
        fprintf(stderr, "Error: CPU count must be a positive integer\n");
        return 0;
    }
    if (cfg.ksize < 5 || cfg.ksize > 32)
    {
        fprintf(stderr, "Error: Only kmer sizes of 5 to 32 are supported\n");
        return 0;
    }
    if (cfg.coverage > 1)
    {
        fprintf(stderr, "Error: Coverage is the proportion of the sequence covered by high kmers and must be between 0 and 1\n");
        return 0;
    }
    if (cfg.depth <= 1)
    {
        fprintf(stderr, "Error: Depth is the number of times a kmer needs to be found before being flagged as high coverage, it must be above 1\n");
        return 0;
    }
    if (cfg.depth < cfg.cpus * 2)
    {
        fprintf(stderr, "Error: Depth must be at least 2 x number of CPUs (for performance reasons; but this version of the program is written to normalise to 50+\n");
        return 0;
    }
    // since we're processing each thread chunk separately.
    cfg.depth = cfg.depth / cfg.cpus;
    cfg.initial_hash_size = (cfg.memory > 0) ? memoryGB2capacity(cfg.memory) : INITIAL_CAPACITY;
    float memory_per_cpu = capacity2memory(cfg.initial_hash_size);
    printf("Initial hash table size set to %'zu (memory ~ %'0.2f Gb for each of %d threads, ~ %'d Gb total))\n", cfg.initial_hash_size, memory_per_cpu, cfg.cpus, cfg.memory);

    if (cfg.initial_hash_size < 100000)
    {
        fprintf(stderr, "Error: initial kmer table size is too small, it should be set to at least 100000 (or leave empty for default)\n");
        return 0;
    }

    return 1;
}

void process_forward_files(char *first_file, char **additional_files)
{
    // Process the first file
    if (access(first_file, R_OK) == 0)
    {
        cfg.forward_file_count++;
        cfg.forward_files = realloc(cfg.forward_files, cfg.forward_file_count * sizeof(char *));
        if (!cfg.forward_files)
        {
            fprintf(stderr, "Memory allocation failed\n");
            exit(1);
        }
        cfg.forward_files[cfg.forward_file_count - 1] = strdup(first_file);
    }
    else
    {
        fprintf(stderr, "Warning: File '%s' does not exist or is not readable. Skipping.\n", first_file);
    }

    // Process additional files
    while (*additional_files != NULL && additional_files[0][0] != '-')
    {
        if (access(*additional_files, R_OK) == 0)
        {
            cfg.forward_file_count++;
            cfg.forward_files = realloc(cfg.forward_files, cfg.forward_file_count * sizeof(char *));
            if (!cfg.forward_files)
            {
                fprintf(stderr, "Memory allocation failed\n");
                exit(1);
            }
            cfg.forward_files[cfg.forward_file_count - 1] = strdup(*additional_files);
        }
        else
        {
            fprintf(stderr, "Warning: File '%s' does not exist or is not readable. Skipping.\n", *additional_files);
        }
        additional_files++;
        optind++;
    }
}

void process_reverse_files(char *first_file, char **additional_files)
{
    // Process the first file
    if (access(first_file, R_OK) == 0)
    {
        cfg.reverse_file_count++;
        cfg.reverse_files = realloc(cfg.reverse_files, cfg.reverse_file_count * sizeof(char *));
        if (!cfg.reverse_files)
        {
            fprintf(stderr, "Memory allocation failed\n");
            exit(1);
        }
        cfg.reverse_files[cfg.reverse_file_count - 1] = strdup(first_file);
    }
    else
    {
        fprintf(stderr, "Warning: File '%s' does not exist or is not readable. Skipping.\n", first_file);
    }

    // Process additional files
    while (*additional_files != NULL && additional_files[0][0] != '-')
    {
        if (access(*additional_files, R_OK) == 0)
        {
            cfg.reverse_file_count++;
            cfg.reverse_files = realloc(cfg.reverse_files, cfg.reverse_file_count * sizeof(char *));
            if (!cfg.reverse_files)
            {
                fprintf(stderr, "Memory allocation failed\n");
                exit(1);
            }
            cfg.reverse_files[cfg.reverse_file_count - 1] = strdup(*additional_files);
        }
        else
        {
            fprintf(stderr, "Warning: File '%s' does not exist or is not readable. Skipping.\n", *additional_files);
        }
        additional_files++;
        optind++;
    }
}

char *create_output_filename(const char *basename, int k, int norm_depth, int thread)
{
    size_t len = strlen(basename) + 30;
    char *output_filename = malloc(len);
    if (!output_filename)
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    snprintf(output_filename, len, "%s.k%d_norm%d_thread%d", basename, k, norm_depth, thread);
    return output_filename;
}

////////////////////////////////////////////////////////////////
////////////////HASH TABLE FUNCTIONS ///////////////////////////
////////////////////////////////////////////////////////////////

bool is_hash_table_empty(const hash_table_t *ht)
{
    for (size_t i = 0; i < ht->capacity; i++)
    {
        if (ht->entries[i].hash != 0)
        {
            return false;
        }
    }
    return true;
}

void init_hash_table(hash_table_t *ht)
{
    ht->entries = calloc(cfg.initial_hash_size, sizeof(kmer_t));
    if (!ht->entries)
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    ht->used = 0;
    ht->capacity = cfg.initial_hash_size;
}

hash_table_t *copy_hash_table(const hash_table_t *source)
{
    hash_table_t *copy = malloc(sizeof(hash_table_t));
    if (!copy)
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    copy->capacity = source->capacity;
    copy->used = source->used;
    copy->entries = calloc(copy->capacity, sizeof(kmer_t));
    if (!copy->entries)
    {
        fprintf(stderr, "Memory allocation failed\n");
        free(copy);
        exit(EXIT_FAILURE);
    }
    memcpy(copy->entries, source->entries, copy->capacity * sizeof(kmer_t));
    return copy;
}

int store_kmer(hash_table_t *hash_table, uint64_t hash, int thread_id)
{
    // this happens locally on a thread specific hash_table so no requirement to do any locking

    if (hash_table->used >= hash_table->capacity * TABLE_LOAD_FACTOR)
        expand_local_hash_table(hash_table, 0, thread_id);

    size_t index = hash % hash_table->capacity;

    // if entry is null and available:
    if (hash_table->entries[index].hash == 0)
    {
        if (cfg.debug > 3)
        {
            char kmer_str[cfg.ksize + 1];
            decode_kmer_plain(hash, cfg.ksize, kmer_str);
            printf("new kmer %s, hash %zu at index %zu. Existing entry is %zu and hash capacity is %zu and used size %zu\n", kmer_str, hash, index, hash_table->entries[index].hash, hash_table->capacity, hash_table->used);
        }
        hash_table->entries[index].hash = hash;
        hash_table->entries[index].count = 1;
        // new kmer inserted, so increase size used.
        hash_table->used++;
        return index;
    }
    else if (hash_table->entries[index].hash == hash)
    {
        hash_table->entries[index].count++;
    }
    else // collision
    {

        size_t original_index = index;
        int collisions = 0;

        while (hash_table->entries[index].hash != 0 && hash_table->entries[index].hash != hash)
        {
            collisions++;
            if (cfg.debug > 2)
                printf("Thread %d: hash %'zu collision consecutive number %d, index: %'zu -> %'zu, capacity %'zu\n", thread_id, hash, collisions, original_index, index, hash_table->capacity);

            if ((collisions / hash_table->capacity) > (hash_table->capacity * 0.5))
            {
                if (cfg.debug)
                    printf("Thread %d: Collisions more than 10%% of table capacity, expanding table...\n", thread_id);
                expand_local_hash_table(hash_table, 0, thread_id);
            }

            index = (index + 1) % hash_table->capacity;

            hash_table->entries[index].count++;
        }
    }

    if (cfg.debug > 3)
        printf("DEBUG: Kmer hash: %lu, Count: %d\n", hash, hash_table->entries[index].count);
    return index;
}

size_t expand_local_hash_table(hash_table_t *ht, size_t new_capacity, int thread_id)
{
    if (!new_capacity || new_capacity == 0)
        new_capacity = ht->capacity + ht->capacity * 0.5;

    if (ht->capacity >= new_capacity)
    {
        return ht->capacity;
    }

    size_t new_hash_size = 0; // not necessary but Just In Case (TM)

    if (cfg.debug)
        printf("Thread %d: Local hash table expansion triggered, from %'zu to %'zu\n", thread_id, ht->capacity, new_capacity);

    kmer_t *new_entries = calloc(new_capacity, sizeof(kmer_t));
    if (!new_entries)
    {
        fprintf(stderr, "Error: Thread %d: Memory allocation failed to expand local hash table, from %'zu to %'zu\n", thread_id, ht->capacity, new_capacity);
        exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < ht->capacity; i++)
    {
        if (ht->entries[i].hash != 0)
        {
            // kmer is stored at an index based on its capacity so when that increases, indexes change
            size_t new_index = ht->entries[i].hash % new_capacity;
            while (new_entries[new_index].hash != 0)
            {
                new_index = (new_index + 1) % new_capacity;
            }
            new_entries[new_index] = ht->entries[i];
            new_hash_size++;
        }
    }

    if (ht->entries)
        free(ht->entries);

    ht->entries = new_entries;
    ht->capacity = new_capacity;
    ht->used = new_hash_size;

    if (ht->used >= ht->capacity * 0.90)
    {
        fprintf(stderr, "Warning: Thread %d: Local hash table is still over 90%% full after expansion (%'zu)\n", thread_id, ht->used);
    }

    if (cfg.debug)
        printf("Thread %d: Local hash table expansion completed successfully, using %'zu (mem: %'.2f) of %'zu new capacity (mem: %'.2f)\n", thread_id, ht->used, capacity2memory(ht->used), ht->capacity, capacity2memory(ht->capacity));

    return ht->capacity;
}

////////////////////////////////////////////////////////////////
//////////////////HASHING AND ENCODING /////////////////////////
////////////////////////////////////////////////////////////////

/////////////////

uint64_t murmurhash3(const char *key, int len, uint32_t seed)
{
    const uint64_t c1 = 0x87c37b91114253d5;
    const uint64_t c2 = 0x4cf5ad432745937f;
    uint64_t h1 = seed;
    const uint64_t *blocks = (const uint64_t *)(key);
    int nblocks = len / 8;

    for (int i = 0; i < nblocks; i++)
    {
        uint64_t k1 = blocks[i];
        k1 *= c1;
        k1 = (k1 << 31) | (k1 >> (64 - 31));
        k1 *= c2;

        h1 ^= k1;
        h1 = (h1 << 27) | (h1 >> (64 - 27));
        h1 = h1 * 5 + 0x52dce729;
    }

    const uint8_t *tail = (const uint8_t *)(key + nblocks * 8);
    uint64_t k1 = 0;

    switch (len & 7)
    {
    case 7:
        k1 ^= ((uint64_t)tail[6]) << 48;
    case 6:
        k1 ^= ((uint64_t)tail[5]) << 40;
    case 5:
        k1 ^= ((uint64_t)tail[4]) << 32;
    case 4:
        k1 ^= ((uint64_t)tail[3]) << 24;
    case 3:
        k1 ^= ((uint64_t)tail[2]) << 16;
    case 2:
        k1 ^= ((uint64_t)tail[1]) << 8;
    case 1:
        k1 ^= ((uint64_t)tail[0]);
        k1 *= c1;
        k1 = (k1 << 31) | (k1 >> (64 - 31));
        k1 *= c2;
        h1 ^= k1;
    }

    h1 ^= len;
    h1 ^= h1 >> 33;
    h1 *= 0xff51afd7ed558ccd;
    h1 ^= h1 >> 33;
    h1 *= 0xc4ceb9fe1a85ec53;
    h1 ^= h1 >> 33;

    return h1;
}

static inline uint64_t encode_kmer_murmur(const char *seq, int k)
{
    // even if the sequence is larger, we only store up to k
    char kmer[k];
    for (int i = 0; i < k; i++)
    {
        kmer[i] = base_map[(uint8_t)seq[i]];
    }

    uint64_t hash1 = murmurhash3(kmer, k, 0);
    // printf("kmer %s hash1 %u hash2 %u\n", seq, hash1, hash2);
    return hash1;
}

/////////////////

static inline uint64_t encode_kmer_plain(const char *kmer, int k)
{
    uint64_t encoded = 0;
    for (int i = 0; i < k; i++)
    {
        encoded = (encoded << 2) | base_map[(uint8_t)kmer[i]];
    }
    return encoded;
}

char decode_kmer_plain(uint64_t encoded, int k, char *kmer)
{
    for (int i = k - 1; i >= 0; i--)
    {
        kmer[i] = rev_base_map[encoded & 0x03];
        encoded >>= 2;
    }
    kmer[k] = '\0';
}

/////////////////

static inline uint64_t mix_bits(uint64_t hash)
{
    hash ^= (hash >> 33);
    hash *= 0xff51afd7ed558ccd;
    hash ^= (hash >> 33);
    hash *= 0xc4ceb9fe1a85ec53;
    hash ^= (hash >> 33);
    return hash;
}

static inline uint64_t kmer_hash_fnv(const char *seq, int k)
{
    uint64_t hash = 0xcbf29ce484222325;
    for (int i = 0; i < k; i++)
    {
        hash ^= base_map[(uint8_t)seq[i]];
        hash *= 0x100000001b3;
    }
    return mix_bits(hash);
}

/////////////////

////////////////////////////////////////////////////////////////
////////////////// DNA MANIPULATION ////////////////////////////
////////////////////////////////////////////////////////////////

bool valid_dna(const char *sequence)
{
    const char *valid_chars = "ATCG";

    for (int i = 0; sequence[i] != '\0'; i++)
    {
        // assume upper case
        if (strchr(valid_chars, sequence[i]) == NULL)
        {
            return false;
        }
    }

    return true;
}

void reverse_complement(const char *seq, char *rev_comp, int k)
{
    // assume no Ns anymore.
    static const char complement[256] = {
        ['A'] = 'T', ['T'] = 'A', ['C'] = 'G', ['G'] = 'C'
        //, ['a'] = 't', ['t'] = 'a', ['c'] = 'g', ['g'] = 'c'
    };

    for (int i = 0; i < k; i++)
        rev_comp[k - 1 - i] = complement[(unsigned char)seq[i]];

    // null term
    rev_comp[k] = '\0';
}

const char *get_canonical_kmer(const char *kmer, int k)
{
    static char rev_comp[MAX_K + 1];
    reverse_complement(kmer, rev_comp, k);
    return (strcmp(kmer, rev_comp) < 0) ? kmer : rev_comp;
}

////////////////////////////////////////////////////////////////
///////////////// DATA PROCESSING //////////////////////////////
////////////////////////////////////////////////////////////////

void process_sequence(const char *seq, hash_table_t *hash_table, int K, int NORM_DEPTH, int *seq_high_count_kmers, int *total_seq_kmers, int thread_id)
{
    *seq_high_count_kmers = 0;
    *total_seq_kmers = 0;
    replacestr(seq, "N", "A");
    int seq_len = strlen(seq);

    if (seq_len >= K && valid_dna(seq))
    {
        for (int i = 0; i <= seq_len - K; i++)
        {
            char kmer[K + 1];
            strncpy(kmer, seq + i, K);
            kmer[K] = '\0';

            // printf("sequence is %s\n", kmer);
            uint64_t hash = 0;

            if (cfg.canonical == 1)
            {
                const char *canonical_kmer = get_canonical_kmer(kmer, K);
                // hash = encode_kmer_murmur(canonical_kmer, K);
                hash = encode_kmer_plain(canonical_kmer, K);
            }
            else
            {
                // hash = encode_kmer_murmur(kmer, K);
                hash = encode_kmer_plain(kmer, K);
            }

            if (hash == 0)
                continue;

            // uint8_t first_base = base_map[(uint8_t)seq[0]];
            // hash_table_t *used_hash_table = &hash_table[first_base];

            (*total_seq_kmers)++;
            int index = store_kmer(hash_table, hash, thread_id);

            if (hash_table->entries[index].count >= NORM_DEPTH)
            {
                (*seq_high_count_kmers)++;
            }
        }
    }
}

void *process_thread_chunk(void *arg)
{
    thread_data_t *data = (thread_data_t *)arg;
    char *forward_ptr = data->forward_ptr;
    char *reverse_ptr = data->reverse_ptr;
    size_t forward_size = data->forward_size;
    size_t reverse_size = data->reverse_size;
    hash_table_t *local_hash_table = data->hash_table;
    int thread_id = data->thread_id;
    FILE *output_forward = data->thread_output_forward;
    FILE *output_reverse = data->thread_output_reverse;

    int lines_to_read = strcmp(cfg.filetype, "fa") == 0 ? 2 : 4;
    char forward_record[lines_to_read][MAX_SEQ_LENGTH], reverse_record[lines_to_read][MAX_SEQ_LENGTH];

    // thread stats
    size_t total_kmers = data->hash_table->used;
    int processed_count = 0, printed_count = 0, skipped_count = 0, prev_printed_count = data->printed_count, prev_skipped_count = data->skipped_count;
    size_t prev_total_kmers = total_kmers, sequences_last_processed = data->processed_count;
    double prev_rate = 0;

    if (cfg.verbose)
    {
        printf("Thread %d started; processing %d lines per record\n", thread_id, lines_to_read);
    }

    size_t sequences_processed = 0;
    time_t current_time = time(NULL);

    // report every 60 seconds

    while (forward_ptr < data->forward_ptr + forward_size && reverse_ptr < data->reverse_ptr + reverse_size)
    {

        bool valid_record = true;
        int seq_high_count_kmers_forward = 0, total_seq_kmers_forward = 0;
        int seq_high_count_kmers_reverse = 0, total_seq_kmers_reverse = 0;
        int seq_index = 1; // line with dna sequence, always 1 in this implementation.

        // read pair
        for (int i = 0; i < lines_to_read; i++)
        {
            forward_ptr = read_line(forward_ptr, forward_record[i], MAX_SEQ_LENGTH);
            reverse_ptr = read_line(reverse_ptr, reverse_record[i], MAX_SEQ_LENGTH);

            if (!forward_ptr || !reverse_ptr)
            {
                valid_record = false;
                break;
            }
        }
        if (!valid_record)
            break;

        if (cfg.debug > 2)
        {
            printf("DEBUG: Thread %d - Processing sequence pair %'zu\n", thread_id, data->processed_count);
            printf("DEBUG: Forward sequence: %.50s...\n", forward_record[seq_index]);
            printf("DEBUG: Reverse sequence: %.50s...\n", reverse_record[seq_index]);
        }

        int forward_seq_len = strlen(forward_record[seq_index]);
        int reverse_seq_len = strlen(reverse_record[seq_index]);
        if (forward_seq_len < cfg.ksize || reverse_seq_len < cfg.ksize)
        {
            data->skipped_count++;
            if (cfg.debug)
            {
                printf("Sequence pair %'zu skipped due to short length (F:%d;R:%d)\n", data->processed_count, forward_seq_len, reverse_seq_len);
            }
            continue; // Skip to the next record
        }
        process_sequence(forward_record[seq_index], local_hash_table, cfg.ksize, cfg.depth, &seq_high_count_kmers_forward, &total_seq_kmers_forward, thread_id);
        process_sequence(reverse_record[seq_index], local_hash_table, cfg.ksize, cfg.depth, &seq_high_count_kmers_reverse, &total_seq_kmers_reverse, thread_id);

        int seq_high_count_kmers = seq_high_count_kmers_forward + seq_high_count_kmers_reverse;
        int total_seq_kmers = total_seq_kmers_forward + total_seq_kmers_reverse;
        float high_count_ratio = total_seq_kmers > 0 ? (float)seq_high_count_kmers / total_seq_kmers : 0;

        data->processed_count++;

        bool is_printed = false;
        if (high_count_ratio <= cfg.coverage)
        {
            for (int i = 0; i < lines_to_read; i++)
            {
                fprintf(output_forward, "%s\n", forward_record[i]);
                fprintf(output_reverse, "%s\n", reverse_record[i]);
            }
            data->printed_count++;
            is_printed = true;
        }
        else
        {
            data->skipped_count++;
            is_printed = false;
        }

        // debug report
        if (cfg.debug > 1)
        {

            if (is_printed == true)
            {
                printf("Thread %d - Sequence pair %'zu PRINTED: ", thread_id, data->processed_count);
            }
            else
            {
                printf("Thread %d - Sequence pair %'zu SKIPPED: ", thread_id, data->processed_count);
            }
            printf("High (%d) count kmers: F:%d;R:%d;B:%d, Total kmers: F:%d;R:%d;B:%d, High (%d) count ratio: %.2f\n",
                   cfg.depth,
                   seq_high_count_kmers_forward,
                   seq_high_count_kmers_reverse, seq_high_count_kmers,
                   total_seq_kmers_forward, total_seq_kmers_reverse, total_seq_kmers,
                   cfg.depth,
                   high_count_ratio);
        }
        sequences_processed++;

        if (difftime(current_time, data->last_report_time) >= REPORTING_INTERVAL)
        {
            if (cfg.debug > 1)
                printf("reporting after %d seconds\n", REPORTING_INTERVAL);
            current_time = time(NULL);
            double elapsed_time = difftime(current_time, data->last_report_time);
            sequences_last_processed = data->processed_count - data->last_report_count;
            double rate = sequences_last_processed / elapsed_time;
            total_kmers = local_hash_table->used;

            if (prev_total_kmers > 0 || prev_printed_count > 0 || prev_skipped_count > 0)
            {
                float printed_improvement = (prev_printed_count == 0) ? 0 : (float)(data->printed_count - prev_printed_count) / prev_printed_count;
                float skipped_improvement = (prev_skipped_count == 0) ? 0 : (float)(data->skipped_count - prev_skipped_count) / prev_skipped_count;
                float prev_rate_improvement = (prev_rate == 0) ? 0 : (float)(rate - prev_rate) / prev_rate;
                float kmer_improvement = (prev_total_kmers == 0) ? 0 : (float)(total_kmers - prev_total_kmers) / prev_total_kmers;
                printf("Thread %d - Processing rate: %'.0f (%+.2f%%) sequences per second, processed %'zu pairs, printed: %'zu (%+.2f%%), skipped: %'zu (%+.2f%%), Total kmers across all sequences: %'zu (%+.2f%%)\n",
                       thread_id, rate, prev_rate_improvement * 100,
                       data->processed_count,
                       data->printed_count,
                       printed_improvement * 100,
                       data->skipped_count,
                       skipped_improvement * 100,
                       total_kmers, kmer_improvement * 100);
            }

            prev_total_kmers = total_kmers;
            prev_printed_count = data->printed_count;
            prev_skipped_count = data->skipped_count;
            prev_rate = rate;

            data->last_report_time = current_time;
            data->last_report_count = data->processed_count;
        }
    }

    if (cfg.verbose)
        printf("Thread %d: completed processing file\n", thread_id);

    current_time = time(NULL);
    double elapsed_time = difftime(current_time, data->last_report_time);
    sequences_processed = data->processed_count - data->last_report_count;
    double rate = sequences_processed / elapsed_time;
    total_kmers = local_hash_table->used;
    float printed_improvement = (prev_printed_count == 0) ? 0 : (float)(data->printed_count - prev_printed_count) / prev_printed_count;
    float skipped_improvement = (prev_skipped_count == 0) ? 0 : (float)(data->skipped_count - prev_skipped_count) / prev_skipped_count;
    float prev_rate_improvement = (prev_rate == 0) ? 0 : (float)(rate - prev_rate) / prev_rate;
    float kmer_improvement = (prev_total_kmers == 0) ? 0 : (float)(total_kmers - prev_total_kmers) / prev_total_kmers;
    printf("Thread %d - Processing rate: %'.0f (%+.2f%%) sequences per second, processed %'zu pairs, printed: %'zu (%+.2f%%), skipped: %'zu (%+.2f%%), Total kmers across all sequences: %'zu (%+.2f%%)\n",
           thread_id, rate, prev_rate_improvement * 100,
           data->processed_count,
           data->printed_count,
           printed_improvement * 100,
           data->skipped_count,
           skipped_improvement * 100,
           total_kmers, kmer_improvement * 100);

    prev_total_kmers = total_kmers;
    prev_printed_count = data->printed_count;
    prev_skipped_count = data->skipped_count;
    prev_rate = rate;

    data->last_report_time = current_time;
    data->last_report_count = data->processed_count;

    return NULL;
}

int multithreaded_process_files(thread_data_t *thread_data, mmap_file_t *forward_mmap, mmap_file_t *reverse_mmap)
{
    char stop_char = strcmp(cfg.filetype, "fa") == 0 ? '>' : '@';

    pthread_t threads[cfg.cpus];

    size_t total_file_size = forward_mmap->size;
    size_t approx_chunk_size = total_file_size / cfg.cpus;
    char *forward_end = forward_mmap->data + total_file_size;
    char *reverse_end = reverse_mmap->data + reverse_mmap->size;

    char *forward_start = forward_mmap->data;
    char *reverse_start = reverse_mmap->data;

    size_t processed_count;
    size_t printed_count;
    size_t skipped_count;
    time_t last_report_time;
    int last_report_count;
    int thread_sync_interval;

    for (int i = 0; i < cfg.cpus; i++)
    {

        if (cfg.debug > 1)
            printf("Starting thread %d\n", i);



        char *forward_chunk_end = (i == cfg.cpus - 1) ? forward_end : find_next_record_start(forward_start + approx_chunk_size, forward_end, stop_char);
        char *reverse_chunk_end = (i == cfg.cpus - 1) ? reverse_end : find_next_record_start(reverse_start + approx_chunk_size, reverse_end, stop_char);

        thread_data[i].forward_ptr = forward_start;
        thread_data[i].reverse_ptr = reverse_start;
        thread_data[i].forward_size = forward_chunk_end - forward_start;
        thread_data[i].reverse_size = reverse_chunk_end - reverse_start;

        if (!thread_data[i].thread_output_forward || !thread_data[i].thread_output_reverse)
        {
            fprintf(stderr, "Error opening thread-specific output files for thread %d, %s and %s\n", i, thread_data[i].thread_output_forward_filename, thread_data[i].thread_output_reverse_filename);
            exit(EXIT_FAILURE);
        }

        // printf("forward/reverse size %zu/%zu chunk %s/%s\n", thread_data[i].forward_size, thread_data[i].reverse_size, forward_chunk_end, reverse_chunk_end);

        if (pthread_create(&threads[i], NULL, process_thread_chunk, &thread_data[i]) != 0)
        {
            fprintf(stderr, "Error creating thread %d\n", i);
            for (int j = 0; j < i; j++)
            {
                pthread_cancel(threads[j]);
                pthread_join(threads[j], NULL);
            }
            return 1;
        }

        forward_start = forward_chunk_end;
        reverse_start = reverse_chunk_end;
        sleep(1);
    }

    // close
    for (int i = 0; i < cfg.cpus; i++)
    {
        if (pthread_join(threads[i], NULL) != 0)
        {
            fprintf(stderr, "Error joining thread %d\n", i);
            return 1;
        }
    }

    // Generate final report
    size_t total_processed = 0, total_printed = 0, total_skipped = 0;
    for (int i = 0; i < cfg.cpus; i++)
    {
        total_processed += thread_data[i].processed_count;
        total_printed += thread_data[i].printed_count;
        total_skipped += thread_data[i].skipped_count;
    }
    reporting.total_processed = total_processed;
    reporting.total_printed = total_printed;
    reporting.total_skipped = total_skipped;
    reporting.files_processed++;

    printf("File's statistics: Processed %'zu, Printed %'zu, Skipped %'zu\n",
           total_processed, total_printed, total_skipped);

    return 0;
}

////////////////////////////////////////////////////////////////
/////////////// MAIN ///////////////////////////////////////////
////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
    setlocale(LC_ALL, "");
    memset(&reporting, 0, sizeof(struct reporting_t));
    reporting.total_processed = 0;
    reporting.total_printed = 0;
    reporting.total_skipped = 0;
    reporting.total_kmers = 0;
    reporting.files_processed = 0;

    memset(&cfg, 0, sizeof(struct config_t));
    if (!parse_arguments(argc, argv))
    {
        printf("Issue parsing options\n");
        print_usage(argv[0]);
        return 1;
    }

    // for (int i=0;i<NUM_PARTITIONS;i++){
    //     init_hash_table(&global_hash_tables[i]);
    // }

    // initialise the output files for each thread (otherwise we'd need to open as append)
    // every thread prints to its own file

    thread_data_t thread_data[cfg.cpus];
    for (int thread_number = 0; thread_number < cfg.cpus; thread_number++)
    {
        thread_data[thread_number].thread_output_forward_filename = NULL;
        thread_data[thread_number].thread_output_reverse_filename = NULL;
        thread_data[thread_number].thread_output_forward_filename = create_output_filename("output_forward", cfg.ksize, cfg.depth, thread_number);
        thread_data[thread_number].thread_output_reverse_filename = create_output_filename("output_reverse", cfg.ksize, cfg.depth, thread_number);

        thread_data[thread_number].thread_output_forward = fopen(thread_data[thread_number].thread_output_forward_filename, "w");
        if (!thread_data[thread_number].thread_output_forward)
        {
            fprintf(stderr, "Error opening file to write: %s\n", thread_data[thread_number].thread_output_forward_filename);
            exit(EXIT_FAILURE);
        }
        thread_data[thread_number].thread_output_reverse = fopen(thread_data[thread_number].thread_output_reverse_filename, "w");
        if (!thread_data[thread_number].thread_output_reverse)
        {
            fprintf(stderr, "Error opening file to write: %s\n", thread_data[thread_number].thread_output_reverse_filename);
            exit(EXIT_FAILURE);
        }

        thread_data[thread_number].thread_sync_interval = SYNC_INTERVAL;
        thread_data[thread_number].thread_id = thread_number;
        thread_data[thread_number].processed_count = 0;
        thread_data[thread_number].printed_count = 0;
        thread_data[thread_number].skipped_count = 0;
        thread_data[thread_number].total_kmers = 0;
        thread_data[thread_number].last_report_time = time(NULL);
        thread_data[thread_number].last_report_count = 0;

        thread_data[thread_number].forward_ptr = NULL;
        thread_data[thread_number].reverse_ptr = NULL;

        thread_data[thread_number].hash_table = malloc(sizeof(hash_table_t));
        if (!thread_data[thread_number].hash_table)
        {
            fprintf(stderr, "Thread %d: Memory allocation failed for hash table\n", thread_number);
            exit(EXIT_FAILURE);
        }
        init_hash_table(thread_data[thread_number].hash_table);

        // thread_data[i].hash_table = &global_hash_table; // pointer,
        // for (int m = 0; i < NUM_PARTITIONS; m++)
        // {
        //     thread_data[i].hash_tables[m] = NULL;
        //     thread_data[i].hash_tables[m] = copy_hash_table(&global_hash_tables[m]);
        // }
        // printf("MAIN2: Is the global table empty? %s\n", is_hash_table_empty(&global_hash_table) ? "true" : "false");
        // printf("MAIN3: Is the local table of thread %d empty? %s\n", i, is_hash_table_empty(thread_data[i].hash_table) ? "true" : "false");
    }
    time_t start_time = time(NULL);

    for (int file_index = 0; file_index < cfg.forward_file_count; file_index++)
    {
        printf("Processing file pair %d of %d: %s and %s\n", file_index + 1, cfg.forward_file_count, cfg.forward_files[file_index], cfg.reverse_files[file_index]);
        mmap_file_t forward_mmap = mmap_file(cfg.forward_files[file_index]);
        mmap_file_t reverse_mmap = mmap_file(cfg.reverse_files[file_index]);

        if (forward_mmap.data == NULL || reverse_mmap.data == NULL)
        {
            fprintf(stderr, "Error memory mapping input files\n");
            goto cleanup;
        }

        // passing local copy as was having issues with pointers.
        if (multithreaded_process_files(thread_data, &forward_mmap, &reverse_mmap) != 0)
        {
            fprintf(stderr, "Error processing files\n");
            munmap_file(&forward_mmap);
            munmap_file(&reverse_mmap);
            goto cleanup;
        }
        munmap_file(&forward_mmap);
        munmap_file(&reverse_mmap);
    }
    for (int thread_number = 0; thread_number < cfg.cpus; thread_number++)
    {
        fclose(thread_data[thread_number].thread_output_forward);
        fclose(thread_data[thread_number].thread_output_reverse);
        free(thread_data[thread_number].thread_output_forward_filename);
        free(thread_data[thread_number].thread_output_reverse_filename);
        free(thread_data[thread_number].hash_table->entries);
        free(thread_data[thread_number].hash_table);
    }

    printf("\n--- Final Report ---\n");
    printf("Processed Records: %'zu\n", reporting.total_processed);
    printf("Printed Records: %'zu\n", reporting.total_printed);
    printf("Skipped Records: %'zu\n", reporting.total_skipped);
    printf("Total kmers across all sequences: %'zu\n", reporting.total_kmers);

cleanup:
    // Clean up
    for (int thread_number = 0; thread_number < cfg.forward_file_count; thread_number++)
    {
        free(cfg.forward_files[thread_number]);
        free(cfg.reverse_files[thread_number]);
    }
    free(cfg.forward_files);
    free(cfg.reverse_files);

    time_t end_time = time(NULL);
    double total_runtime = difftime(end_time, start_time);
    printf("Program completed, cleaning up.\n");
    printf("Total runtime: %.2f seconds\n", total_runtime);
    // if (total_processed > 0)
    // {
    //     double processing_rate = total_processed / total_runtime;
    //     printf("Overall processing rate: %.2f sequences per second\n", processing_rate);
    // }
    // else
    // {
    //     printf("No data processed\n");
    // }
    return 0;
}
