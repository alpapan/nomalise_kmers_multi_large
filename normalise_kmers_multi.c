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

// TODO: provide an option to store kmers canonically. it would have tye memory needed too

// fix memory.leak when new file is rewd

// #define INITIAL_CAPACITY 500000003 // equiv genome/transcriptome size; should be a prime number; remember we do not convert to canonical kmers (as this was create for rnaseq)
#define INITIAL_CAPACITY 150000001 // good starting value for a transcriptome, esp if stranded
#define MAX_SEQ_LENGTH 400
#define CHUNK_SIZE 1000000
#define MAX_THREADS 256
#define SYNC_INTERVAL 100000 // how many seqs/often to sync the kmer tables of each thread, this inits the thread-specific sync interval which increases with more data.

typedef struct hash_table_t hash_table_t;

typedef struct
{
    char *forward_ptr;
    char *reverse_ptr;
    size_t forward_size;
    size_t reverse_size;
    hash_table_t *hash_table;
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
// pthread_mutex_t output_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t sync_table_mutex = PTHREAD_MUTEX_INITIALIZER;

typedef struct
{
    uint64_t hash;
    int count;
} kmer_t;

struct hash_table_t
{
    kmer_t *entries;
    size_t size;
    size_t capacity;
};

struct reporting
{
    size_t total_processed;
    size_t total_printed;
    size_t total_skipped;
    size_t total_kmers;
};

struct config
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
void expand_global_hash_table(hash_table_t *ht, size_t new_capacity, int thread_id);
void expand_local_hash_table(hash_table_t *ht, size_t new_capacity, int thread_id);
static inline uint64_t kmer_hash(const char *seq, int k);
char *create_output_filename(const char *basename, int k, int norm_depth, int thread);
int multithreaded_process_files(thread_data_t *thread_data, mmap_file_t *forward_mmap, mmap_file_t *reverse_mmap, hash_table_t *hash_table);
int find_or_insert_kmer(hash_table_t *hash_table, uint64_t hash, int thread_id);
void process_sequence(const char *seq, hash_table_t *hash_table, int K, int NORM_DEPTH, int *seq_high_count_kmers, int *total_seq_kmers, int thread_id);
void print_usage(char *program_name);
int parse_arguments(int argc, char *argv[]);
mmap_file_t mmap_file(const char *filename);
void munmap_file(mmap_file_t *mf);
void synchronise_hash_tables(hash_table_t *local_ht, hash_table_t *global_ht, int thread_id);
void *process_thread_chunk(void *arg);
bool valid_dna(const char *sequence);
char *read_line(char *ptr, char *buffer, int max_length);
static void replacestr(const char *line, const char *search, const char *replace);
char *find_next_record_start(char *ptr, char *end_ptr, char stop_char);

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
    fprintf(stderr, "Usage: %s\t--forward file1 [file2+] --reverse file1 [file2+] [--ksize (int; def. 25)] [--depth|-d (int; def. 100)] \n\t\t\t\t[--coverage|g (float 0-1; def. 0.9)] [--verbose] [--filetype|-t (fq|fa; def. fq)] [--cpu|-p (int; def 4)] [--debug|-b]\n", program_name);
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
    cfg.filetype = "fq";
    cfg.cpus = 4;
    cfg.help = 0;
    cfg.forward_file_count = 0;
    cfg.reverse_file_count = 0;
    cfg.forward_files = NULL;
    cfg.reverse_files = NULL;
    cfg.ksize = 25;
    cfg.depth = 100;

    static struct option long_options[] = {
        {"forward", required_argument, 0, 'f'},
        {"reverse", required_argument, 0, 'r'},
        {"ksize", required_argument, 0, 'k'},
        {"depth", required_argument, 0, 'd'},
        {"coverage", required_argument, 0, 'g'},
        {"filetype", required_argument, 0, 't'},
        {"cpu", required_argument, 0, 'p'},
        {"debug", required_argument, 0, 'b'},
        {"verbose", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}};

    int opt;
    int option_index = 0;

    while ((opt = getopt_long(argc, argv, "f:r:k:d:g:t:p:b:vh", long_options, &option_index)) != -1)
    {
        switch (opt)
        {
        case 'b':
            cfg.debug = atoi(optarg);
            break;
        case 'h':
            cfg.help = 1;
            return 0;
        case 'p':
            cfg.cpus = atoi(optarg);
            if (cfg.cpus <= 0)
            {
                fprintf(stderr, "Error: CPU count must be a positive integer\n");
                return 0;
            }
            break;
        case 'f':
            process_forward_files(optarg, &argv[optind]);

            break;
        case 'r':
            process_reverse_files(optarg, &argv[optind]);

            break;
        case 'k':
            cfg.ksize = atoi(optarg);
            if (cfg.ksize < 5 || cfg.ksize > 32)
            {
                fprintf(stderr, "Error: Only kmer sizes of 5 to 32 are supported\n");
                return 0;
            }
            break;
        case 'd':
            cfg.depth = atoi(optarg);
            if (cfg.depth <= 1)
            {
                fprintf(stderr, "Error: Depth is the number of times a kmer needs to be found before being flagged as high coverage, it must be above 1\n");
                return 0;
            }
            break;
        case 'g':
            cfg.coverage = atof(optarg);
            if (cfg.coverage > 1)
            {
                fprintf(stderr, "Error: Coverage is the proportion of the sequence covered by high kmers and must be between 0 and 1\n");
                return 0;
            }
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

    return 1;
}

//////////////////////////////


void init_hash_table(hash_table_t *ht)
{
    pthread_mutex_lock(&hash_table_mutex);
    ht->entries = calloc(INITIAL_CAPACITY, sizeof(kmer_t));
    if (!ht->entries)
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    ht->size = 0;
    ht->capacity = INITIAL_CAPACITY;
    pthread_mutex_unlock(&hash_table_mutex);
}

void expand_global_hash_table(hash_table_t *ht, size_t new_capacity, int thread_id)
{
    if (!new_capacity || new_capacity == 0)
        new_capacity = ht->capacity + ht->capacity * 0.5; // increase by 50%

    if (ht->capacity >= new_capacity)
        return;

    pthread_mutex_lock(&hash_table_mutex);
    if (cfg.debug)
        printf("Thread %d: Global hash table expansion triggered, from %'zu to %'zu\n", thread_id, ht->capacity, new_capacity);

    kmer_t *new_entries = calloc(new_capacity, sizeof(kmer_t));
    if (!new_entries)
    {
        fprintf(stderr, "Error: Thread %d: Memory allocation failed to expand global hash table, from %'zu to %'zu\n", thread_id, ht->capacity, new_capacity);
        exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < ht->capacity; i++)
    {
        if (ht->entries[i].hash != 0)
        {
            size_t new_index = ht->entries[i].hash % new_capacity;
            while (new_entries[new_index].hash != 0)
            {
                new_index = (new_index + 1) % new_capacity;
            }
            new_entries[new_index] = ht->entries[i];
        }
    }

    free(ht->entries);
    ht->entries = new_entries;
    ht->capacity = new_capacity;

    if (ht->size >= ht->capacity * 0.90)
    {
        fprintf(stderr, "Warning: Thread %d: Global hash table is still over 90%% full after expansion (%'zu)\n", thread_id, ht->size);
    }

    if (cfg.debug)
        printf("Thread %d: Global hash table expansion completed successfully, using %'zu of %'zu new_capacity\n", thread_id, ht->size, ht->capacity);

    pthread_mutex_unlock(&hash_table_mutex);
}

void expand_local_hash_table(hash_table_t *ht, size_t new_capacity, int thread_id)
{
    if (!new_capacity || new_capacity == 0)
        new_capacity = ht->capacity + ht->capacity * 0.5;

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
            size_t new_index = ht->entries[i].hash % new_capacity;
            while (new_entries[new_index].hash != 0)
            {
                new_index = (new_index + 1) % new_capacity;
            }
            new_entries[new_index] = ht->entries[i];
        }
    }

    free(ht->entries);
    ht->entries = new_entries;
    ht->capacity = new_capacity;

    if (ht->size >= ht->capacity * 0.90)
    {
        fprintf(stderr, "Warning: Thread %d: Local hash table is still over 90%% full after expansion (%'zu)\n", thread_id, ht->size);
    }

    if (cfg.debug)
        printf("Thread %d: Local hash table expansion completed successfully, using %'zu of %'zu new_capacity\n", thread_id, ht->size, ht->capacity);
}

static inline uint64_t kmer_hash(const char *seq, int k)
{
    uint64_t hash = 0;
    for (int i = 0; i < k; i++)
    {
        hash = (hash << 2) | base_map[(uint8_t)seq[i]];
    }
    return hash;
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

int find_or_insert_kmer(hash_table_t *hash_table, uint64_t hash, int thread_id)
{
    // this happens locally on a thread specific hash_table so no requirement to do any locking

    if (hash_table->size >= hash_table->capacity * 0.9)
    {
        expand_local_hash_table(hash_table, 0, thread_id);
    }

    size_t index = hash % hash_table->capacity;
    size_t original_index = index;

    while (hash_table->entries[index].hash != 0 && hash_table->entries[index].hash != hash)
    {
        index = (index + 1) % hash_table->capacity;
        if (index == original_index)
        {
            if (cfg.debug)
                printf("hash table is full at index %'zu. This shouldn't happen really but will try to increase it.\n", index);
            expand_local_hash_table(hash_table, 0, thread_id);
            return find_or_insert_kmer(hash_table, hash, thread_id);
        }
    }

    if (hash_table->entries[index].hash == 0)
    {
        hash_table->entries[index].hash = hash;
        hash_table->entries[index].count = 0;
        hash_table->size++;
    }

    hash_table->entries[index].count++;

    //   printf("DEBUG: Kmer hash: %lu, Count: %d\n", hash, hash_table->entries[index].count);
    return index;
}

void process_sequence(const char *seq, hash_table_t *hash_table, int K, int NORM_DEPTH, int *seq_high_count_kmers, int *total_seq_kmers, int thread_id)
{
    *seq_high_count_kmers = 0;
    *total_seq_kmers = 0;
    replacestr(seq, "N", "A");
    int seq_len = strlen(seq);

    if (seq_len > K && valid_dna(seq))

        for (int i = 0; i <= seq_len - K; i++)
        {
            uint64_t hash = kmer_hash(seq + i, K);
            if (hash == 0)
                continue;

            (*total_seq_kmers)++;
            int index = find_or_insert_kmer(hash_table, hash, thread_id);

            if (hash_table->entries[index].count >= NORM_DEPTH)
            {
                (*seq_high_count_kmers)++;
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
    hash_table_t *hash_table = data->hash_table;
    int thread_id = data->thread_id;
    FILE *output_forward = data->thread_output_forward;
    FILE *output_reverse = data->thread_output_reverse;

    int lines_to_read = strcmp(cfg.filetype, "fa") == 0 ? 2 : 4;
    char forward_record[lines_to_read][MAX_SEQ_LENGTH], reverse_record[lines_to_read][MAX_SEQ_LENGTH];

    // thread stats
    int processed_count = 0, printed_count = 0, skipped_count = 0, prev_printed_count = 0, prev_skipped_count = 0;
    size_t total_kmers = 0, prev_total_kmers = 0;
    double prev_rate = 0;

    int sync_frequency_round = 0;

    if (cfg.verbose)
    {
        printf("Thread %d started; processing %d lines per record\n", thread_id, lines_to_read);
    }

    hash_table_t local_hash_table;
    init_hash_table(&local_hash_table);

    size_t sequences_processed = 0;

    data->last_report_time = time(NULL);
    data->last_report_count = 0;
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
        process_sequence(forward_record[seq_index], &local_hash_table, cfg.ksize, cfg.depth, &seq_high_count_kmers_forward, &total_seq_kmers_forward, thread_id);
        process_sequence(reverse_record[seq_index], &local_hash_table, cfg.ksize, cfg.depth, &seq_high_count_kmers_reverse, &total_seq_kmers_reverse, thread_id);

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
            printf("High (%d) count kmers: F:%d;R:%d;B:%d, Total kmers: F:%d;R:%d;B:%d, High (%d) count ratio: %.4f\n",
                   cfg.depth,
                   seq_high_count_kmers_forward,
                   seq_high_count_kmers_reverse, seq_high_count_kmers,
                   total_seq_kmers_forward, total_seq_kmers_reverse, total_seq_kmers,
                   cfg.depth,
                   high_count_ratio);
        }
        sequences_processed++;

        // report every 60 seconds
        time_t current_time = time(NULL);
        if (difftime(current_time, data->last_report_time) >= 60)
        {
            double elapsed_time = difftime(current_time, data->last_report_time);
            size_t sequences_processed = data->processed_count - data->last_report_count;
            double rate = sequences_processed / elapsed_time;
            total_kmers = local_hash_table.size;
            float printed_improvement = (prev_printed_count == 0) ? 0 : (float)(data->printed_count - prev_printed_count) / prev_printed_count;
            float skipped_improvement = (prev_skipped_count == 0) ? 0 : (float)(data->skipped_count - prev_skipped_count) / prev_skipped_count;
            float prev_rate_improvement = (prev_rate == 0) ? 0 : (float)(rate - prev_rate) / prev_rate;
            float kmer_improvement = (prev_total_kmers == 0) ? 0 : (float)(total_kmers - prev_total_kmers) / prev_total_kmers;

            // once we hit a million skipped and skipped > printed then double the data->thread_sync_interval but only once
            // TODO: perhaps this should be a function of kmer_improvement, at some point we are not seeing many new kmers.
            if (sync_frequency_round == 0 && kmer_improvement < 0.10)
            {
                sync_frequency_round++;
                data->thread_sync_interval = SYNC_INTERVAL * 2;
                if (cfg.debug)
                    printf("Thread %d: Increasing current sync frequency by 10x to %d\n", thread_id, data->thread_sync_interval);
            }
            else if (sync_frequency_round == 1 && data->skipped_count > 1e6 && data->skipped_count > data->printed_count)
            {
                sync_frequency_round++;
                data->thread_sync_interval = SYNC_INTERVAL * 2;
                if (cfg.debug)
                    printf("Thread %d: Increasing current sync frequency by 5x to %d\n", thread_id, data->thread_sync_interval);
            }
            else if (sync_frequency_round == 2 && data->skipped_count > 1e8 && data->skipped_count > data->printed_count && kmer_improvement < 0.05)
            {
                sync_frequency_round++;
                data->thread_sync_interval = SYNC_INTERVAL * 10;
                if (cfg.debug)
                    printf("Thread %d: Increasing current sync frequency by 2x to %d\n", thread_id, data->thread_sync_interval);
            }

            printf("Thread %d - Processing rate: %.0f (%+.2f%%) sequences per second, processed %'zu pairs, printed: %'zu (%+.2f%%), skipped: %'zu (%+.2f%%), Total kmers across all sequences: %'zu (%+.2f%%)\n",
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
        }

        if (sequences_processed % data->thread_sync_interval == 0)
        {
            synchronise_hash_tables(&local_hash_table, hash_table, thread_id);
        }
    }
    synchronise_hash_tables(&local_hash_table, hash_table, thread_id);

    if (cfg.verbose)
        printf("Thread %d: completed processing file\n", thread_id);

    return NULL;
}

int multithreaded_process_files(thread_data_t *thread_data, mmap_file_t *forward_mmap, mmap_file_t *reverse_mmap, hash_table_t *global_hash_table)
{
    char stop_char = strcmp(cfg.filetype, "fa") == 0 ? '>' : '@';

    pthread_t threads[cfg.cpus];

    size_t total_size = forward_mmap->size;
    size_t approx_chunk_size = total_size / cfg.cpus;
    char *forward_end = forward_mmap->data + total_size;
    char *reverse_end = reverse_mmap->data + reverse_mmap->size;

    char *forward_start = forward_mmap->data;
    char *reverse_start = reverse_mmap->data;

    char *forward_ptr;
    char *reverse_ptr;
    size_t forward_size;
    size_t reverse_size;
    hash_table_t *hash_table;
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


	// pointers must be freed?

         thread_data[i].forward_ptr = NULL;
         thread_data[i].reverse_ptr = NULL;
//keep         thread_data[i].hash_table = NULL;

    }

    // Generate final report
    size_t total_processed = 0, total_printed = 0, total_skipped = 0;
    for (int i = 0; i < cfg.cpus; i++)
    {
        total_processed += thread_data[i].processed_count;
        total_printed += thread_data[i].printed_count;
        total_skipped += thread_data[i].skipped_count;
    }
    // reported_data->total_processed = total_processed;
    // reported_data->total_printed = total_printed;
    // reported_data->total_skipped = total_skipped;
    // reported_data->total_kmers = total_kmers;

    printf("File's statistics: Processed %'zu, Printed %'zu, Skipped %'zu, Total unique kmers: %'zu\n",
           total_processed, total_printed, total_skipped, global_hash_table->size);

    return 0;
}

void synchronise_hash_tables(hash_table_t *local_ht, hash_table_t *global_ht, int thread_id)
{
    // only one thread can sync at a time
    // pthread_mutex_lock(&sync_table_mutex);

    // check if hash tables are the same size, otherwise make them so.
    if (local_ht->capacity > global_ht->capacity)
    {
        if (cfg.debug)
            printf("Thread %d: Increasing global table capacity %'zu to meet local %'zu\n", thread_id, global_ht->capacity, local_ht->capacity);
        expand_global_hash_table(global_ht, local_ht->capacity, thread_id);
    }
    if (global_ht->capacity > local_ht->capacity)
    {
        if (cfg.debug)
            printf("Thread %d: Increasing local table capacity %'zu to meet global %'zu\n", thread_id, local_ht->capacity, global_ht->capacity);
        expand_local_hash_table(local_ht, global_ht->capacity, thread_id);
    }

    if (cfg.debug)
        printf("Thread %d: Synchronising local hash table (%'zu) with master (%'zu)\n", thread_id, local_ht->capacity, global_ht->capacity);

    if (global_ht->capacity != local_ht->capacity)
    {
        if (cfg.debug)
            printf("Thread %d: Global hash table capacity %'zu is not the same as local table %'zu!\n", thread_id, global_ht->capacity, local_ht->capacity);
        exit(EXIT_FAILURE);
    }

    if (cfg.debug)
        printf("Thread %d: Locking global hash table\n", thread_id);

    pthread_mutex_lock(&hash_table_mutex);
    // Merge local hash table into global hash table
    for (size_t i = 0; i < local_ht->capacity; i++)
    {
        if (local_ht->entries[i].hash != 0)
        {
            size_t index = local_ht->entries[i].hash % global_ht->capacity;
            while (global_ht->entries[index].hash != 0 &&
                   global_ht->entries[index].hash != local_ht->entries[i].hash)
            {
                index = (index + 1) % global_ht->capacity;
            }
            if (global_ht->entries[index].hash == 0)
            {
                global_ht->entries[index] = local_ht->entries[i];
                __atomic_add_fetch(&global_ht->size, 1, __ATOMIC_SEQ_CST);
            }
            else
            {
                __atomic_add_fetch(&global_ht->entries[index].count,
                                   local_ht->entries[i].count, __ATOMIC_SEQ_CST);
            }
        }
    }

    // Copy global hash table back to local hash table
    free(local_ht->entries);
    local_ht->entries = calloc(global_ht->capacity, sizeof(kmer_t));
    if (!local_ht->entries)
    {
        fprintf(stderr, "Thread %d: Memory allocation failed in synchronize_hash_tables\n", thread_id);
        exit(EXIT_FAILURE);
    }
    memcpy(local_ht->entries, global_ht->entries, global_ht->capacity * sizeof(kmer_t));
    local_ht->size = global_ht->size;
    local_ht->capacity = global_ht->capacity;

    if (cfg.debug)
        printf("Thread %d: Unlocked global hash table\n", thread_id);
    pthread_mutex_unlock(&hash_table_mutex);
    // pthread_mutex_unlock(&sync_table_mutex);

    if (cfg.debug)
        printf("Thread %d: Hash table sync complete\n", thread_id);
}

int main(int argc, char *argv[])
{
    setlocale(LC_ALL, "");
    memset(&cfg, 0, sizeof(struct config));

    if (!parse_arguments(argc, argv))
    {
        printf("Issue parsing options\n");
        print_usage(argv[0]);
        return 1;
    }

    // Initialize mutexes
    if (pthread_mutex_init(&hash_table_mutex, NULL) != 0 ||
        // pthread_mutex_init(&output_mutex, NULL) != 0 ||
        pthread_mutex_init(&sync_table_mutex, NULL) != 0)
    {
        fprintf(stderr, "Error initializing mutexes\n");
        return 1;
    }

    // char *output_forward_filename = create_output_filename("output_forward", cfg.ksize, cfg.depth);
    // char *output_reverse_filename = create_output_filename("output_reverse", cfg.ksize, cfg.depth);

    // FILE *output_forward = fopen(output_forward_filename, "w");
    // FILE *output_reverse = fopen(output_reverse_filename, "w");
    // if (!output_forward || !output_reverse)
    // {
    //     fprintf(stderr, "Error opening output files\n");
    //     goto cleanup;
    // }

    hash_table_t global_hash_table;
    init_hash_table(&global_hash_table);

    // initialise the output files for each thread (otherwise we'd need to open as append)
    // every thread prints to its own file
    thread_data_t thread_data[cfg.cpus];
    for (int i = 0; i < cfg.cpus; i++)
    {
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
        thread_data[i].hash_table = &global_hash_table; // copy global hash table into local thread

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
        if (multithreaded_process_files(thread_data, &forward_mmap, &reverse_mmap, &global_hash_table) != 0)
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
        //free(thread_data[i].thread_output_reverse);
        //free(thread_data[i].thread_output_reverse);
 }
    printf("\n--- Final Report ---\n");
    // printf("Processed Records: %'zu\n", reporting->total_processed);
    // printf("Printed Records: %'zu\n", reporting->total_printed);
    // printf("Skipped Records: %'zu\n", reporting->total_skipped);
    // printf("Total kmers across all sequences: %'zu\n", reporting->total_kmers);
    printf("Final hash table size: %'zu\n", global_hash_table.size);

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
    // pthread_mutex_destroy(&output_mutex);
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

