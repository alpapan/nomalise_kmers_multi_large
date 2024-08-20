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

// use plain hash
// organise
// remove redumdant code
// support trinity output
// ensure records are read 2/4 lines at a time
// support gz bz2 input will be hard because of record boundaries
// support gz bz2 output will be easier

// #define INITIAL_CAPACITY 500000003 // 7.5gb per cpu equiv genome/transcriptome size; should be a prime number; remember we do not convert to canonical kmers (as this was create for rnaseq)
#define INITIAL_CAPACITY 150000001 // 2.3gb good starting value for a transcriptome, esp if stranded
#define MAX_LINE_LENGTH 1024
#define CHUNK_SIZE 1000000
#define MAX_THREADS 256
// #define SYNC_INTERVAL 100000 // how many seqs/often to sync the kmer tables of each thread, this inits the thread-specific sync interval which increases with more data.
#define SYNC_INTERVAL 1000000
#define TABLE_LOAD_FACTOR 0.6
#define MAX_K 32
#define REPORTING_INTERVAL 60 // seconds

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
hash_table_t global_hash_table;

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

pthread_mutex_t hash_table_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t sync_table_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t reporting_data = PTHREAD_MUTEX_INITIALIZER;

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
    char *informat;
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

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

static const uint8_t base_map[256] = {
    ['A'] = 0x00, ['C'] = 0x01, ['G'] = 0x02, ['T'] = 0x03};

// Function prototypes
void init_hash_table(hash_table_t *ht);
size_t expand_global_hash_table(size_t new_capacity, int thread_id, bool nolock);
size_t expand_local_hash_table(hash_table_t *ht, size_t new_capacity, int thread_id);
static inline uint64_t encode_kmer_murmur(const char *seq, int k);
char *create_output_filename(const char *basename, int k, int norm_depth, int thread);
int multithreaded_process_files(thread_data_t *thread_data, mmap_file_t *forward_mmap, mmap_file_t *reverse_mmap);
size_t store_kmer(hash_table_t *hash_table, uint64_t hash, int thread_id);
void process_sequence(const char *seq, hash_table_t *hash_table, int K, int NORM_DEPTH, int *seq_high_count_kmers, int *total_seq_kmers, int thread_id);
void print_usage(char *program_name);
int parse_arguments(int argc, char *argv[]);
mmap_file_t mmap_file(const char *filename);
void munmap_file(mmap_file_t *mf);
void synchronise_hash_tables(hash_table_t *local_ht, int thread_id);
void *process_thread_chunk(void *arg);
bool valid_dna(const char *sequence);
char *read_line(char *ptr, char *buffer, int max_length);
static void replacestr(const char *line, const char *search, const char *replace);
char *find_next_record_start(char *ptr, char *end_ptr, char stop_char);
hash_table_t *copy_hash_table(const hash_table_t *source);
float capacity2memory(size_t capacity);
float memoryGB2capacity(size_t memory);
uint64_t find_hash_offset(uint64_t hash, size_t capacity);

////////////////
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
//////////////////////////////
void print_usage(char *program_name)
{
    fprintf(stderr, "Usage: %s\t--forward file1 [file2+] --reverse file1 [file2+] [--ksize (int; def. 25)] [--depth|-d (int; def. 100)]\n\
    \t\t\t[--coverage|g (float 0-1; def. 0.9)] [--verbose] [--filetype|-t (fq|fa; def. fq)] [--cpu|-p (int; def 4)] [--debug|-b]\n\
    [--canonical|c] [--memory|m (int; def. 150000001)] \n",
            program_name);
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

int parse_arguments(int argc, char *argv[])
{
    cfg.coverage = 0.9;
    cfg.verbose = 0;
    cfg.debug = 0;
    cfg.informat = "fq";
    cfg.cpus = 4;
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
            cfg.informat = optarg;
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

//////////////////////////////
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

void init_hash_table(hash_table_t *ht)
{
    pthread_mutex_lock(&hash_table_mutex);
    ht->entries = calloc(cfg.initial_hash_size, sizeof(kmer_t));
    if (!ht->entries)
    {
        fprintf(stderr, "Memory allocation failed\n");
        pthread_mutex_unlock(&hash_table_mutex);
        exit(EXIT_FAILURE);
    }
    ht->used = 0;
    ht->capacity = cfg.initial_hash_size;
    pthread_mutex_unlock(&hash_table_mutex);
}

hash_table_t *copy_hash_table(const hash_table_t *source)
{
    pthread_mutex_lock(&hash_table_mutex);
    hash_table_t *copy = malloc(sizeof(hash_table_t));
    if (!copy)
    {
        fprintf(stderr, "Memory allocation failed\n");
        pthread_mutex_unlock(&hash_table_mutex);
        exit(EXIT_FAILURE);
    }
    copy->capacity = source->capacity;
    copy->used = source->used;
    copy->entries = calloc(copy->capacity, sizeof(kmer_t));
    if (!copy->entries)
    {
        fprintf(stderr, "Memory allocation failed\n");
        free(copy);
        pthread_mutex_unlock(&hash_table_mutex);
        exit(EXIT_FAILURE);
    }
    memcpy(copy->entries, source->entries, copy->capacity * sizeof(kmer_t));
    pthread_mutex_unlock(&hash_table_mutex);
    return copy;
}
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
size_t expand_global_hash_table(size_t new_capacity, int thread_id, bool nolock)
{
    // no need to worry about global hash table if 1 cpu.
    if (cfg.cpus == 1)
        return 0;

    if (!new_capacity || new_capacity == 0)
        new_capacity = global_hash_table.capacity + global_hash_table.capacity * 0.5; // increase by 50%

    if (global_hash_table.capacity >= new_capacity)
    {
        return global_hash_table.capacity;
    }

    if (cfg.debug)
    {
        printf("Thread %d: Master Hash table expansion triggered, from %'zu to %'zu\n", thread_id, global_hash_table.capacity, new_capacity);
    }

    kmer_t *new_entries = calloc(new_capacity, sizeof(kmer_t));
    if (!new_entries)
    {
        fprintf(stderr, "Error: Thread %d: Memory allocation failed to expand Master Hash table, from %'zu to %'zu\n", thread_id, global_hash_table.capacity, new_capacity);
        exit(EXIT_FAILURE);
    }

    size_t new_global_hash_size = 0; // not necessary but Just In Case (TM)

    if (nolock == false)
        pthread_mutex_lock(&hash_table_mutex);

    for (size_t i = 0; i < global_hash_table.capacity; i++)
    {
        if (global_hash_table.entries[i].hash != 0)
        {
            // kmer is stored at an index based on its capacity so when that increases, indexes change
            size_t new_index = global_hash_table.entries[i].hash % new_capacity;
            uint64_t offset = find_hash_offset(global_hash_table.entries[i].hash, new_capacity);

            while (new_entries[new_index].hash != 0)
            {
                new_index = (new_index + offset) % new_capacity;
            }
            new_entries[new_index] = global_hash_table.entries[i];

            new_global_hash_size++;
        }
    }

    if (global_hash_table.entries)
        free(global_hash_table.entries);

    global_hash_table.entries = new_entries;
    global_hash_table.capacity = new_capacity;
    global_hash_table.used = new_global_hash_size;

    if (new_global_hash_size >= new_capacity * 0.90)
    {
        fprintf(stderr, "Warning: Thread %d: Master Hash table is still over 90%% full after expansion (%'zu)\n", thread_id, new_global_hash_size);
    }

    if (cfg.debug)
        printf("Thread %d: Master Hash table expansion completed successfully, using %'zu (mem: %'.2f) of %'zu new capacity (mem: %'.2f)\n", thread_id, new_global_hash_size, capacity2memory(new_global_hash_size), new_capacity, capacity2memory(new_capacity));

    if (nolock == false)
        pthread_mutex_unlock(&hash_table_mutex);

    return global_hash_table.capacity;
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
            uint64_t offset = find_hash_offset(ht->entries[i].hash, new_capacity);
            while (new_entries[new_index].hash != 0)
            {
                new_index = (new_index + offset) % new_capacity;
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

///////////////////////////

uint64_t find_hash_offset(uint64_t hash, size_t capacity)
{

    if (!capacity || capacity == 0)
       capacity = global_hash_table.capacity;

    return 1 + (hash % (capacity - 1));
}

/////////////////////////

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
    char kmer[k];
    for (int i = 0; i < k; i++)
    {
        kmer[i] = base_map[(uint8_t)seq[i]];
    }

    uint64_t hash1 = murmurhash3(kmer, k, 0);
    // printf("kmer %s hash1 %u hash2 %u\n", seq, hash1, hash2);
    return hash1;
}

//////////////////////////

static inline uint64_t encode_kmer_plain(const char *seq, int k)
{
    uint64_t hash = 0;
    for (int i = 0; i < k; i++)
    {
        hash = (hash << 2) | base_map[(uint8_t)seq[i]];
    }
    return hash;
}

///////////////////

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

//////////////////////

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

size_t store_kmer(hash_table_t *hash_table, uint64_t hash, int thread_id)
{
    // this happens locally on a thread specific hash_table so no requirement to do any locking

    if (hash_table->used >= hash_table->capacity * TABLE_LOAD_FACTOR)
    {
        expand_local_hash_table(hash_table, 0, thread_id);
    }

    size_t index = hash % hash_table->capacity;
    if (index < 0)
    {
        fprintf(stderr, "ERROR1: This shouldnt have happened, index %zu hash %zu capacity %zu\n", index, hash, hash_table->capacity);
        exit(EXIT_FAILURE);
    }

    size_t original_index = index;
    uint64_t offset = find_hash_offset(hash, hash_table->capacity);

    while (hash_table->entries[index].hash != 0 && hash_table->entries[index].hash != hash)
    {
        index = (index + offset) % hash_table->capacity;
        if (index < 0)
        {
            fprintf(stderr, "ERROR2: This shouldnt have happened, index %zu hash %zu capacity %zu\n", index, hash, hash_table->capacity);
            exit(EXIT_FAILURE);
        }

        if (index == 0)
        {
            if (cfg.debug)
                printf("hash table is full at index %'zu. This shouldn't happen really but will try to increase it.\n", index);
            expand_local_hash_table(hash_table, 0, thread_id);
            return store_kmer(hash_table, hash, thread_id);
        }
    }

    // new kmer inserted, so increase size.
    if (hash_table->entries[index].hash == 0)
    {
        hash_table->entries[index].hash = hash;
        hash_table->entries[index].count = 0;
        hash_table->used++;
    }

    hash_table->entries[index].count++;

    //   printf("DEBUG: Kmer hash: %lu, Count: %d\n", hash, hash_table->entries[index].count);
    return index;
}

void reverse_complement(const char *seq, char *rev_comp, int k)
{
    static const char complement[256] = {
        ['A'] = 'T', ['T'] = 'A', ['C'] = 'G', ['G'] = 'C'};

    for (int i = 0; i < k; i++)
        rev_comp[k - 1 - i] = complement[(unsigned char)seq[i]];

    // null termi
    rev_comp[k] = '\0';
}

const char *get_canonical_kmer(const char *kmer, int k)
{
    static char rev_comp[MAX_K + 1];
    reverse_complement(kmer, rev_comp, k);
    return (strcmp(kmer, rev_comp) < 0) ? kmer : rev_comp;
}

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
                hash = encode_kmer_murmur(canonical_kmer, K);
            }
            else
            {
                hash = encode_kmer_murmur(kmer, K);
            }

            if (hash == 0)
                continue;

            // uint8_t first_base = base_map[(uint8_t)seq[0]];
            // hash_table_t *used_hash_table = &hash_table[first_base];

            (*total_seq_kmers)++;
            size_t index = store_kmer(hash_table, hash, thread_id);
            if (index < 0)
            {
                fprintf(stderr, "ERROR1: This shouldnt have happened, index %zu hash %zu capacity %zu\n", index, hash, hash_table->capacity);
                exit(EXIT_FAILURE);
            }

            if (hash_table->entries[index].count >= NORM_DEPTH)
            {
                (*seq_high_count_kmers)++;
            }
        }
    }
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

// faster and checks for both stop char and number of lines (thanks gpt!)
char *adjust_chunk_to_record_boundary(char *start_ptr, char *end_ptr, size_t approx_chunk_size, char stop_char)
{
    char *chunk_end_ptr = start_ptr + approx_chunk_size;
    int lines_to_read = strcmp(cfg.informat, "fa") == 0 ? 2 : 4;

    if (chunk_end_ptr > end_ptr)
    {
        chunk_end_ptr = end_ptr;
    }

    size_t lines_count = 0;
    char *ptr = start_ptr;

    // Iterate over the data until the end of the chunk or end of data
    while (ptr < chunk_end_ptr)
    {
        // Use a direct pointer increment to find the next newline
        while (ptr < chunk_end_ptr && *ptr != '\n')
        {
            ptr++;
        }
        if (ptr < chunk_end_ptr)
        {
            lines_count++;
            ptr++; // Move past the newline
            // Update chunk_end_ptr only when a complete record is found
            if (lines_count % lines_to_read == 0)
            {
                chunk_end_ptr = ptr;
            }
        }
    }

    return chunk_end_ptr;
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

    int lines_to_read = strcmp(cfg.informat, "fa") == 0 ? 2 : 4;
    char forward_record[lines_to_read][MAX_LINE_LENGTH], reverse_record[lines_to_read][MAX_LINE_LENGTH];

    // thread stats
    size_t total_kmers = data->hash_table->used;
    int processed_count = 0, printed_count = 0, skipped_count = 0, prev_printed_count = data->printed_count, prev_skipped_count = data->skipped_count;
    size_t prev_total_kmers = total_kmers, sequences_last_processed = data->processed_count;
    double prev_rate = 0;

    int sync_frequency_round = 0;

    if (cfg.verbose)
    {
        printf("Thread %d started; processing %d lines per record\n", thread_id, lines_to_read);
    }

    size_t sequences_processed = 0;
    time_t current_time = time(NULL);

    while (forward_ptr < data->forward_ptr + forward_size && reverse_ptr < data->reverse_ptr + reverse_size)
    {

        bool valid_record = true;
        int seq_high_count_kmers_forward = 0, total_seq_kmers_forward = 0;
        int seq_high_count_kmers_reverse = 0, total_seq_kmers_reverse = 0;
        int seq_index = 1; // line with dna sequence, always 1 in this implementation.

        // read pair
        for (int i = 0; i < lines_to_read; i++)
        {
            forward_ptr = read_line(forward_ptr, forward_record[i], MAX_LINE_LENGTH);
            reverse_ptr = read_line(reverse_ptr, reverse_record[i], MAX_LINE_LENGTH);

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
            // pthread_mutex_lock(&output_mutex); // no longer needed as we have 1 file per thread
            for (int i = 0; i < lines_to_read; i++)
            {
                fprintf(output_forward, "%s\n", forward_record[i]);
                fprintf(output_reverse, "%s\n", reverse_record[i]);
            }
            // pthread_mutex_unlock(&output_mutex);
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
            printf("High (%d) count kmers: F:%d;R:%d;B:%d, Total unique kmers: F:%d;R:%d;B:%d, High (%d) count ratio: %.4f\n",
                   cfg.depth,
                   seq_high_count_kmers_forward,
                   seq_high_count_kmers_reverse, seq_high_count_kmers,
                   total_seq_kmers_forward, total_seq_kmers_reverse, total_seq_kmers,
                   cfg.depth,
                   high_count_ratio);
        }
        sequences_processed++;

        // report every 60 seconds
        current_time = time(NULL);
        if (difftime(current_time, data->last_report_time) >= 60)
        {
            double elapsed_time = difftime(current_time, data->last_report_time);
            size_t sequences_processed = data->processed_count - data->last_report_count;
            double rate = sequences_processed / elapsed_time;
            total_kmers = local_hash_table->used;

            float printed_improvement = (prev_printed_count == 0) ? 0 : (float)(data->printed_count - prev_printed_count) / prev_printed_count;
            float skipped_improvement = (prev_skipped_count == 0) ? 0 : (float)(data->skipped_count - prev_skipped_count) / prev_skipped_count;
            float prev_rate_improvement = (prev_rate == 0) ? 0 : (float)(rate - prev_rate) / prev_rate;
            float kmer_improvement = (prev_total_kmers == 0) ? 0 : (float)(total_kmers - prev_total_kmers) / prev_total_kmers;
            printf("Thread %d - Processing rate: %'.0f (%+.2f%%) sequences per second, processed %'zu pairs, printed: %'zu (%+.2f%%), skipped: %'zu (%+.2f%%), Total unique kmers across all sequences: %'zu (%+.2f%%)\n",
                   thread_id, rate, prev_rate_improvement * 100,
                   data->processed_count,
                   data->printed_count,
                   printed_improvement * 100,
                   data->skipped_count,
                   skipped_improvement * 100,
                   total_kmers, kmer_improvement * 100);

            // once we hit a million skipped and skipped > printed then double the data->thread_sync_interval but only once
            // TODO: perhaps this should be a function of kmer_improvement, at some point we are not seeing many new kmers.
            if (sync_frequency_round == 0 && kmer_improvement < 0.10)
            {
                sync_frequency_round++;
                data->thread_sync_interval = SYNC_INTERVAL * 2;
                if (cfg.debug)
                    printf("Thread %d: Increasing current sync frequency by 10x to %'d\n", thread_id, data->thread_sync_interval);
            }
                else if (sync_frequency_round == 1 && data->skipped_count > 1e6 && data->skipped_count > data->printed_count)
                {
                    sync_frequency_round++;
                    data->thread_sync_interval = SYNC_INTERVAL * 2;
                    if (cfg.debug)
                        printf("Thread %d: Increasing current sync frequency by 5x to %'d\n", thread_id, data->thread_sync_interval);
                }
                else if (sync_frequency_round == 2 && data->skipped_count > 1e8 && data->skipped_count > data->printed_count && kmer_improvement < 0.05)
                {
                    sync_frequency_round++;
                    data->thread_sync_interval = SYNC_INTERVAL * 10;
                    if (cfg.debug)
                        printf("Thread %d: Increasing current sync frequency by 2x to %'d\n", thread_id, data->thread_sync_interval);
                }
            prev_total_kmers = total_kmers;
            prev_printed_count = data->printed_count;
            prev_skipped_count = data->skipped_count;
            prev_rate = rate;

            data->last_report_time = current_time;
            data->last_report_count = data->processed_count;
        }

        if (cfg.cpus > 1 && sequences_processed % data->thread_sync_interval == 0)
        {
            synchronise_hash_tables(local_hash_table, thread_id);
        }
    }

    //  if (cfg.cpus > 1)
    synchronise_hash_tables(local_hash_table, thread_id);

    current_time = time(NULL);
    double elapsed_time = difftime(current_time, data->last_report_time);
    sequences_processed = data->processed_count - data->last_report_count;
    double rate = sequences_processed / elapsed_time;
    total_kmers = local_hash_table->used;
    float printed_improvement = (prev_printed_count == 0) ? 0 : (float)(data->printed_count - prev_printed_count) / prev_printed_count;
    float skipped_improvement = (prev_skipped_count == 0) ? 0 : (float)(data->skipped_count - prev_skipped_count) / prev_skipped_count;
    float prev_rate_improvement = (prev_rate == 0) ? 0 : (float)(rate - prev_rate) / prev_rate;
    float kmer_improvement = (prev_total_kmers == 0) ? 0 : (float)(total_kmers - prev_total_kmers) / prev_total_kmers;
    printf("Thread %d - Processing rate: %'.0f (%+.2f%%) sequences per second, processed %'zu pairs, printed: %'zu (%+.2f%%), skipped: %'zu (%+.2f%%), Total unique kmers across all sequences: %'zu (%+.2f%%)\n",
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

    if (cfg.verbose)
        printf("Thread %d: completed processing file\n", thread_id);

    return NULL;
}

int multithreaded_process_files(thread_data_t *thread_data, mmap_file_t *forward_mmap, mmap_file_t *reverse_mmap)
{
    char stop_char = strcmp(cfg.informat, "fa") == 0 ? '>' : '@';

    pthread_t threads[cfg.cpus];

    size_t total_size = forward_mmap->size;
    size_t approx_chunk_size = total_size / cfg.cpus;
    char *forward_end = forward_mmap->data + total_size;
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
        sleep(3);
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
    pthread_mutex_lock(&reporting_data); // we don't have threads anymore but still.
    reporting.total_processed = total_processed;
    reporting.total_printed = total_printed;
    reporting.total_skipped = total_skipped;
    reporting.total_kmers = global_hash_table.used;
    reporting.files_processed++;
    pthread_mutex_unlock(&reporting_data);

    printf("File's statistics: Processed %'zu, Printed %'zu, Skipped %'zu, Total unique kmers: %'zu\n",
           total_processed, total_printed, total_skipped, global_hash_table.used);

    return 0;
}

void synchronise_hash_tables(hash_table_t *local_ht, int thread_id)
{

    // no need to worry about global hash table if 1 cpu.
    if (cfg.cpus == 1)
        return;

    if (cfg.debug > 2)
        printf("Thread %d: Is the global table empty? %s\n", thread_id, is_hash_table_empty(&global_hash_table) ? "true" : "false");

    // check if hash tables are the same size, otherwise make them so.
    // only needed for global table since local table is recreated
    // if (global_hash_table.capacity > local_ht->capacity)
    // {
    //     if (cfg.debug)
    //         printf("Thread %d: Increasing local table capacity %'zu to meet global %'zu\n", thread_id, local_ht->capacity, global_hash_table.capacity);
    //     expand_local_hash_table(local_ht, global_hash_table.capacity, thread_id);
    // }
    if (local_ht->capacity > global_hash_table.capacity)
    {
        if (cfg.debug)
            printf("Thread %d: Increasing global table capacity %'zu to meet local %'zu\n", thread_id, global_hash_table.capacity, local_ht->capacity);

        expand_global_hash_table(local_ht->capacity, thread_id, false);
    }

    if (cfg.debug)
        printf("Thread %d: Synchronising local hash table (%'zu) with master (%'zu)\n", thread_id, local_ht->capacity, global_hash_table.capacity);

    if (global_hash_table.capacity < local_ht->capacity)
    {
        if (cfg.debug)
            printf("Thread %d: Master Hash table capacity %'zu is (still) smaller than the local table %'zu!\n", thread_id, global_hash_table.capacity, local_ht->capacity);

        exit(EXIT_FAILURE);
    }

    if (cfg.debug > 2)
        printf("Thread %d: Locking Master Hash table\n", thread_id);

    pthread_mutex_lock(&hash_table_mutex);

    // Merge (potentially smaller) local hash table into Master Hash table

    size_t new_global_size = global_hash_table.used;
    int collision_count = 0, total_collisions = 0;
    int collision_cutoff = global_hash_table.capacity * 0.5;

    if (cfg.debug > 2)
        printf("Thread %d: up to %d collisions allowed into global table with capacity %zu and used size %zu\n", thread_id, collision_cutoff, global_hash_table.capacity, global_hash_table.used);

    for (size_t i = 0; i < local_ht->capacity; i++)
    {
        // printf("Thread %d: Checking index %zu\n", thread_id, i);
        if (local_ht->entries[i].hash != 0)
        {
            // printf("Thread %d: Kmer found at index %zu\n", thread_id, i);

            // index is hash vakue modulo table.it is going to be stored in
            size_t index = local_ht->entries[i].hash % global_hash_table.capacity;

            // it's empty, set and go to next.
            if (global_hash_table.entries[index].hash == 0)
            {
                global_hash_table.entries[index] = local_ht->entries[i];
                new_global_size++; // we update size later, once.
                continue;
            }

            size_t original_index = index;
            uint64_t offset = find_hash_offset(local_ht->entries[i].hash, 0);
            // printf("Thread %d: Master Hash at index %zu has hash %lu and local is %lu\n", thread_id, index, global_hash_table.entries[index].hash, local_ht->entries[i].hash);

            while (global_hash_table.entries[index].hash != 0 &&
                   global_hash_table.entries[index].hash != local_ht->entries[i].hash)
            {
                if (cfg.debug > 2)
                    if (collision_count == 1)
                        printf("Thread %d: Master Hash is %lu for index %zu, local is %lu, collision %d\n", thread_id, global_hash_table.entries[index].hash, index, local_ht->entries[i].hash, collision_count);

                index = (index + offset) % global_hash_table.capacity;
//                  index = (index + 1) % global_hash_table.capacity;

                collision_count++;
                // thread tables can host different kmers, when they sync the global may run out of capacity.
                if (index == 0 || collision_count != 0 && collision_count > collision_cutoff)
                {
                    if (cfg.debug > 2)
                        printf("Thread %d: Collisions %d are above limit of %d, expanding Master Hash table\n", thread_id, collision_count, collision_cutoff);
                    expand_global_hash_table(0, thread_id, true);
                    total_collisions += collision_count;
                    collision_count = 0;
                    collision_cutoff = global_hash_table.capacity * 0.1;
                }
            }

            if (global_hash_table.entries[index].hash == 0)
            {
                global_hash_table.entries[index] = local_ht->entries[i];
                new_global_size++; // we update size later, once.
            }
            else
            {
                __atomic_add_fetch(&global_hash_table.entries[index].count,
                                   local_ht->entries[i].count, __ATOMIC_SEQ_CST);
            }
        }
    }

    total_collisions += collision_count;

    // set size once
    global_hash_table.used = new_global_size;

    // Copy Master Hash table back to local hash table
    if (local_ht->entries)
        free(local_ht->entries);
    local_ht->entries = calloc(global_hash_table.capacity, sizeof(kmer_t));
    if (!local_ht->entries)
    {
        fprintf(stderr, "Thread %d: Memory allocation failed in synchronize_hash_tables\n", thread_id);
        exit(EXIT_FAILURE);
    }

    memcpy(local_ht->entries, global_hash_table.entries, global_hash_table.capacity * sizeof(kmer_t));
    local_ht->used = global_hash_table.used;
    local_ht->capacity = global_hash_table.capacity;

    if (cfg.debug > 2)
        printf("Thread %d: Unlocked Master Hash table\n", thread_id);
    pthread_mutex_unlock(&hash_table_mutex);
    if (cfg.debug)
        printf("Thread %d: Hash table sync complete, %'d collisions occurred\n", thread_id, total_collisions);
}

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

    // Initialize mutexes
    if (pthread_mutex_init(&hash_table_mutex, NULL) != 0 ||
        pthread_mutex_init(&reporting_data, NULL) != 0 ||
        pthread_mutex_init(&sync_table_mutex, NULL) != 0)
    {
        fprintf(stderr, "Error initializing mutexes\n");
        return 1;
    }

    // for (int i=0;i<NUM_PARTITIONS;i++){
    //     init_hash_table(&global_hash_tables[i]);
    // }
    init_hash_table(&global_hash_table);

    // initialise the output files for each thread (otherwise we'd need to open as append)
    // every thread prints to its own file

    thread_data_t thread_data[cfg.cpus];
    for (int i = 0; i < cfg.cpus; i++)
    {
        thread_data[i].thread_output_forward_filename = NULL;
        thread_data[i].thread_output_reverse_filename = NULL;
        thread_data[i].thread_output_forward_filename = create_output_filename("output_forward", cfg.ksize, cfg.depth, i);
        thread_data[i].thread_output_reverse_filename = create_output_filename("output_reverse", cfg.ksize, cfg.depth, i);
        thread_data[i].thread_output_forward = fopen(thread_data[i].thread_output_forward_filename, "w");
        thread_data[i].thread_output_reverse = fopen(thread_data[i].thread_output_reverse_filename, "w");
        thread_data[i].thread_sync_interval = SYNC_INTERVAL;
        thread_data[i].thread_id = i;
        thread_data[i].processed_count = 0;
        thread_data[i].printed_count = 0;
        thread_data[i].skipped_count = 0;
        thread_data[i].total_kmers = 0;
        thread_data[i].last_report_time = time(NULL);
        thread_data[i].last_report_count = 0;

        thread_data[i].forward_ptr = NULL;
        thread_data[i].reverse_ptr = NULL;

        // copy Master Hash table into local thread
        thread_data[i].hash_table = malloc(sizeof(hash_table_t));
        if (!thread_data[i].hash_table)
        {
            fprintf(stderr, "Thread %d: Memory allocation failed for hash table\n", i);
            exit(EXIT_FAILURE);
        }
        init_hash_table(thread_data[i].hash_table);
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
   for (int i = 0; i < cfg.cpus; i++)
    {
        fclose(thread_data[i].thread_output_forward);
        fclose(thread_data[i].thread_output_reverse);
        free(thread_data[i].thread_output_forward_filename);
        free(thread_data[i].thread_output_reverse_filename);
        free(thread_data[i].hash_table->entries);
        free(thread_data[i].hash_table);
    }

    printf("\n--- Final Report ---\n");
    printf("Processed Records: %'zu\n", reporting.total_processed);
    printf("Printed Records: %'zu\n", reporting.total_printed);
    printf("Skipped Records: %'zu\n", reporting.total_skipped);
    printf("Total unique kmers across all sequences: %'zu\n", reporting.total_kmers);

cleanup:
    // Clean up
    free(global_hash_table.entries);
    for (int i = 0; i < cfg.forward_file_count; i++)
    {
        free(cfg.forward_files[i]);
        free(cfg.reverse_files[i]);
    }
    free(cfg.forward_files);
    free(cfg.reverse_files);

    // Destroy mutexes
    pthread_mutex_destroy(&hash_table_mutex);
    pthread_mutex_destroy(&reporting_data);
    pthread_mutex_destroy(&sync_table_mutex);

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

