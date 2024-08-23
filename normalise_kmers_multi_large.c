#define VERSION 20240823
// NB: written with some AI assistance

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
// DONE to think: a compromise for not using a shared kmer table would be to
// use 1 thread to process the first X (2-3?) million pairs of the first(?) file pair to seed the table
// that wouldn't need to parse the entire file, just process until Xm records are done
// it should take about 3-5 minutes at 10k reads per sec.
// then copy that kmer table to every other thread.
// for simplicity and consistency (and speed by avoiding write i/o), we could
// not print any output during seeding and just process the entire file again
// (it will have to be parsed for line counting read anyway)
// 20240822: this has now been done

// DONE: support single end data
// 20240822: this has now been done

// TODO: add advice of kmer size vs memory, for example a k=15 with 16gb ram per hash will produce no collisions.
// calc 16*4^15/1024^3=16
// but speed increase is not always obvious as lower k requires more kmers calculated
// 20240823: this has now been done

// TODO: initiate a table sync when there is no consequence on wait time:
// when the file is being read to decide where to split it amongst the threads

// TODO: provide advice that input data should be adaptor/quality trimmed
// noting that i haven't tested that (and non-trimmed files will be faster to process
// due to fwd/rev file sizes/record boundaries being identical).

// TODO: provide advice that two rounds of normalisations may be better especially if have lots of files.
// split the work into two sets between two computers and then run it a second time with more stringent --coverage

// TODO: some peculiar behaviour found after double processin
// i used 10 threads and --coverage 0.96 and then 0.98 to process 3.4Tb -> 500 Gb ->
// the second round's thread 0 and 1 found quite a bit of redundancy (0 was 33% and 1 was 15%) but not other threads.
// more concerning is that the number of unique kmers were much lower in thread0/1.
// i have no idea why this could be so since each thread is independent processing a portion of a file:
//
// Thread 0 - Processing rate: 17,950 (+0.00%) sequences/s, processed 16,676,541 pairs, printed: 9,807,759 (+2.81%), skipped: 6,868,782 (+13.34%), Unique kmers (all sequences; this thread): 145,068,177 (+2.06%)
// Thread 1 - Processing rate: 17,042 (+0.00%) sequences/s, processed 16,622,086 pairs, printed: 12,601,054 (+3.37%), skipped: 4,021,032 (+17.93%), Unique kmers (all sequences; this thread): 187,553,106 (+2.31%)
// Thread 2 - Processing rate: 15,679 (+0.00%) sequences/s, processed 16,540,264 pairs, printed: 15,043,236 (+4.45%), skipped: 1,497,028 (+25.07%), Unique kmers (all sequences; this thread): 240,241,727 (+2.10%)
// Thread 3 - Processing rate: 15,229 (+0.00%) sequences/s, processed 16,513,297 pairs, printed: 15,889,814 (+5.36%), skipped: 623,483 (+20.21%), Unique kmers (all sequences; this thread): 274,340,123 (+2.52%)
// Thread 4 - Processing rate: 15,551 (+0.00%) sequences/s, processed 16,532,616 pairs, printed: 15,526,649 (+5.15%), skipped: 1,005,967 (+20.71%), Unique kmers (all sequences; this thread): 274,321,928 (+2.79%)
// Thread 5 - Processing rate: 16,315 (+0.00%) sequences/s, processed 16,578,460 pairs, printed: 16,331,745 (+5.93%), skipped: 246,715 (+35.79%), Unique kmers (all sequences; this thread): 294,211,911 (+2.77%)
// Thread 6 - Processing rate: 16,117 (+0.00%) sequences/s, processed 16,566,560 pairs, printed: 16,290,330 (+6.07%), skipped: 276,230 (+14.73%), Unique kmers (all sequences; this thread): 305,058,106 (+2.96%)
// Thread 7 - Processing rate: 15,301 (+0.00%) sequences/s, processed 16,517,591 pairs, printed: 16,427,108 (+5.84%), skipped: 90,483 (+13.79%), Unique kmers (all sequences; this thread): 324,452,671 (+2.58%)
// Thread 8 - Processing rate: 15,645 (+0.00%) sequences/s, processed 16,538,238 pairs, printed: 16,508,447 (+5.96%), skipped: 29,791 (+46.80%), Unique kmers (all sequences; this thread): 332,330,298 (+2.76%)
// Thread 9 - Processing rate: 16,930 (+0.00%) sequences/s, processed 16,615,328 pairs, printed: 16,580,786 (+6.45%), skipped: 34,542 (+44.09%), Unique kmers (all sequences; this thread): 330,700,470 (+2.99%)
// this only happens for files produced by the later threads... what is different between threads except file coordinates?
// and it's only to do with the output, not input
//
//

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

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

// [ ] - AP 20240821 TODO: remove redundant and debugging code
// [x] - i attempted to integrated a job queue system but then realised it wouldn't actually improve performance;
//        file is still wholly read while converting records to bytes and sending them to the queue.
//        it would only help if the file was not memory mapped as it wouldn't need to store data in memory. but with mmap there is no benefit?
// [ ] - support gz bz2 input will be hard because of mmapped record boundaries, only solution would be to create a temp file?
// [ ] - support streaming gz output will be easier but streaming bz2 not possible (unless compress at the end but then what's the point)?

//////// HELPERS

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

/// registry of allocs (malloc realloc calloc)
// output_filename
// hash_table_t
// hash_table_t->kmer
// forward_pointer // managed by mmap
// reverse_pointer // managed by mmap
// cfg.forward_files
// cfg.reverse_files
// forward_thread_starts
// forward_thread_ends
// reverse_thread_starts
// reverse_thread_ends

// #define INITIAL_CAPACITY 150000001 // 2.3gb good starting value for a transcriptome, esp if stranded;
//  500000003: 7.5gb per cpu equiv genome/transcriptome size; remember we do not convert to canonical kmers (as this was create for rnaseq)
//  we moved away from specifying capacity in code, user specifies starting memory.
#define INITIAL_CAPACITY 67108879 // prime for 1 gb

#define MAX_LINE_LENGTH 1024
#define MAX_OUTPUT_LENGTH 8192
#define CHUNK_SIZE 1e6
#define MAX_THREADS 256
#define TABLE_LOAD_FACTOR 0.8
#define MAX_K 32
#define REPORTING_INTERVAL 60 // seconds
#define SEED_NUMBER 3e6

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
    kmer_t *kmer;
    size_t used;
    size_t capacity;
} hash_table_t;

// if we decide on partitions.
// #define NUM_PARTITIONS 4
// hash_table_t global_hash_tables[NUM_PARTITIONS];

typedef struct
{
    hash_table_t *hash_table;
    // hash_table_t *hash_tables[NUM_PARTITIONS];
    int thread_id;
    size_t processed_count;
    size_t printed_count;
    size_t skipped_count;
    time_t last_report_time;
    int last_report_count;
    FILE *thread_output_forward;
    FILE *thread_output_reverse;
    char *thread_output_forward_filename;
    char *thread_output_reverse_filename;
    char *forward_data;
    char *reverse_data;
    size_t forward_file_start;
    size_t forward_file_end;
    size_t reverse_file_start;
    size_t reverse_file_end;
} thread_data_t;

thread_data_t thread_data[MAX_THREADS];

struct reporting_t
{
    size_t total_processed;
    size_t total_printed;
    size_t total_skipped;
    size_t max_total_kmers;
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
    int depth_per_cpu;
    float coverage;
    int verbose;
    char *informat;
    bool is_input_fastq;
    char *outformat;
    bool is_output_fastq;
    int cpus;
    int help;
    int debug;
    int canonical;
    int memory;
    size_t initial_hash_size;
    int version;
    int singleend;
    int doprint;
} cfg;

typedef struct
{
    char *data;
    size_t size;
} mmap_file_t;

// Function prototypes
size_t to_power(size_t base, size_t exp);
int compare_kmer_count_asc(const void *a, const void *b);
int compare_kmer_count_desc(const void *a, const void *b);
void sort_kmer_table(hash_table_t *ht, bool do_desc, int thread_id);
void print_kmer_table(hash_table_t *ht, char *identifier, bool do_desc, int thread_id);
size_t sizet_increment_even(size_t num);
char *read_line(char *ptr, char *buffer, int max_length);
float capacity2memory(size_t capacity);
size_t memoryGB2capacity(int memory);
mmap_file_t mmap_file(const char *filename);
void munmap_file(mmap_file_t *mf);
static void replacestr(const char *line, const char *search, const char *replace);
////

void print_usage();
int parse_arguments(int argc, char *argv[]);
void process_forward_files(char *first_file, char **additional_files);
void process_reverse_files(char *first_file, char **additional_files);
char *create_output_filename(const char *basename, int k, int norm_depth, int thread, char *suffix);
void fastq_to_fasta(char (*fastq_record)[4][MAX_LINE_LENGTH], char *fasta_output, bool is_forward);

////
void destroy_hash_table(hash_table_t *ht);
void init_hash_table(hash_table_t *ht);
hash_table_t *copy_hash_table(const hash_table_t *source);
size_t store_kmer(hash_table_t *hash_table, uint64_t hash, int thread_id, bool do_init, char *kmer);
size_t expand_local_hash_table(hash_table_t *ht, size_t new_capacity, int thread_id);

////
static inline uint64_t encode_kmer_plain(const char *seq, int k);
char decode_kmer_plain(uint64_t encoded, int k, char *kmer);

////
bool valid_dna(const char *sequence);
void reverse_complement(const char *seq, char *rev_comp, int k);
const char *get_canonical_kmer(const char *kmer, int k);

////
size_t find_thread_exact_end(char *data, size_t start_pos, size_t end_pos);
bool process_sequence_pair(const char *fwd_seq, const char *rev_seq, hash_table_t *hash_table, int *fwd_seq_high_count_kmers,
                           int *fwd_total_seq_kmers, int *rev_seq_high_count_kmers, int *rev_total_seq_kmers, int thread_id, int K, int depth_per_cpu);
void calculate_thread_positions(char *data, size_t total_file_size, int thread_count, size_t **thread_starts, size_t **thread_ends);
void calculate_thread_positions_from_records(char *data, size_t total_file_size, int thread_count, size_t **thread_starts, size_t **thread_ends, size_t records);
size_t count_records_seqfile(char *data, size_t size);
void seed_kmer_hash(char *file, int records_to_seed, hash_table_t *ht, int K, int depth_per_cpu);
void sequence_to_hash(const char *seq, hash_table_t *hash_table, int *seq_high_count_kmers, int *total_seq_kmers, int K, int depth_per_cpu, int thread_id);
void sequence_to_hash_zero(const char *seq, hash_table_t *hash_table, int K, int depth_per_cpu);
bool is_valid_sequence_single(const char *fwd_seq, int K);
bool is_valid_sequence_pair(const char *fwd_seq, const char *rev_seq, int K);
void *process_thread_chunk_paired(void *arg);
int multithreaded_process_files_paired(thread_data_t *thread_data, mmap_file_t *forward_mmap, mmap_file_t *reverse_mmap);
void *process_thread_chunk_single(void *arg);
int multithreaded_process_files_single(thread_data_t *thread_data, mmap_file_t *forward_mmap);

////////////////////////////////////////////////////////////////
////////// HELPER FUNCTIONS ////////////////////////////////////
////////////////////////////////////////////////////////////////
size_t to_power(size_t base, size_t exp)
{
    size_t result = 1;
    while (exp)
    {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }
    return result;
}

int compare_kmer_asc(const void *a, const void *b)
{
    kmer_t *kmerA = (kmer_t *)a;
    kmer_t *kmerB = (kmer_t *)b;
    char kmerA_str[cfg.ksize + 1];
    char kmerB_str[cfg.ksize + 1];
    decode_kmer_plain(kmerA->hash, cfg.ksize, kmerA_str);
    decode_kmer_plain(kmerB->hash, cfg.ksize, kmerB_str);
    return (*(char *)kmerA_str - *(char *)kmerB_str);
}

int compare_kmer_count_asc(const void *a, const void *b)
{
    kmer_t *kmerA = (kmer_t *)a;
    kmer_t *kmerB = (kmer_t *)b;

    return kmerA->count - kmerB->count;
}

int compare_kmer_count_desc(const void *a, const void *b)
{
    kmer_t *kmerA = (kmer_t *)a;
    kmer_t *kmerB = (kmer_t *)b;

    return kmerB->count - kmerA->count;
}

void sort_kmer_table(hash_table_t *ht, bool do_desc, int thread_id)
{

    if (cfg.verbose)
    {
        if (thread_id < 0)
            printf("Seed: sorting hash table\n");
        else
            printf("Thread %d - sorting hash table\n", thread_id);
    }
    // if (do_desc == true)
    //     qsort(ht->kmer, ht->used, sizeof(kmer_t), compare_kmer_count_desc);
    // else
    //     qsort(ht->kmer, ht->used, sizeof(kmer_t), compare_kmer_count_asc);
    // qsort(ht->kmer, ht->capacity, sizeof(kmer_t), compare_kmer_asc);
}

void print_kmer_table(hash_table_t *ht, char *identifier, bool do_desc, int thread_id)
{
    // sort_kmer_table(ht, do_desc, thread_id);

    if (cfg.verbose)
    {
        if (thread_id < 0)
            printf("Seed: printing kmer table\n");
        else
            printf("Thread %d - printing kmer table\n", thread_id);
    }
    char name[24];
    snprintf(name, sizeof(name), "output_kmer%s", identifier);
    char *kmer_file = create_output_filename(name, cfg.ksize, cfg.depth_per_cpu, thread_id, "tsv");
    FILE *kmer_out = fopen(kmer_file, "w");

    for (size_t index = 0; index < ht->capacity; index++)
    {

        if (ht->kmer[index].hash && ht->kmer[index].hash != 0)
        {
            char kmer_str[cfg.ksize + 1];
            decode_kmer_plain(ht->kmer[index].hash, cfg.ksize, kmer_str);

            if (kmer_str)
                fprintf(kmer_out, "%s\t%d\n", kmer_str, ht->kmer[index].count);
        }
    }

    fclose(kmer_out);
    free(kmer_file);
}

size_t sizet_increment_even(size_t num)
{
    if (num % 2 == 0)
        return num + 1;

    return num;
}
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
        ptr++;
    }

    return (*ptr == '\0') ? NULL : ptr;
}

float capacity2memory(size_t capacity)
{
    return (float)capacity * 16 / 1073741824;
}

size_t memoryGB2capacity(int memory)
{
    size_t number = (float)memory * 1073741824;
    float total_mem = number / 16;
    size_t per_cpu = total_mem / cfg.cpus; // not a prime but a bit better than having an even number
    return sizet_increment_even(per_cpu);
}

mmap_file_t mmap_file(const char *filename)
{
    mmap_file_t result = {NULL, 0};
    int fd = open(filename, O_RDONLY);
    if (fd == -1)
    {
        perror("Error opening file");
        return result;
    }

    struct stat sb;
    if (fstat(fd, &sb) == -1)
    {
        perror("Error getting file size");
        close(fd);
        return result;
    }

    result.size = sb.st_size;
    if (cfg.debug > 1)
        printf("File size: %zu\n", result.size);

    result.data = mmap(NULL, result.size, PROT_READ, MAP_PRIVATE, fd, 0);
    if (result.data == MAP_FAILED)
    {
        perror("Error mapping file");
        result.data = NULL;
        result.size = 0;
    }
    else
    {
        if (cfg.debug > 1)
            printf("File mapped successfully\n");

        close(fd);
        return result;
    }
}

void munmap_file(mmap_file_t *mf)
{
    if (mf->data != NULL)
    {
        munmap(mf->data, mf->size);
        mf->data = NULL;
        mf->size = 0;
    }
    if (cfg.debug > 1)
        printf("File unmapped successfully\n");
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

void print_usage()
{
    fprintf(stderr, "Usage:"
                    "\n\n\t\tMandatory:"
                    "\n\t\t* --forward|-f file1 [file2+]\tList of forward (read1) sequence files"
                    "\n\t\t* --reverse|-r file1 [file2+]\tList of reverse (read2) sequence files"
                    "\n\n\t\tOptional:"
                    "\n\t\t[--single|-s] data are single ended, any --forward files not matched with --reverse will be treated as single-end"
                    "\n\t\t[--ksize|-k (integer 5-32; def. 25)]\tNumber of what size of K to use (must be between 5 and 32)"
                    "\n\t\t[--depth|-d (integer; def. 100)]\tNumber determining when a kmer is tagged as high coverage (defaults to 100),"
                    "\n\t\t\t\t\t\t\tmust be above 2xCPU count as each CPU calculates depth independently"
                    "\n\t\t[--coverage|-g (float 0-1; def. 0.9)]\tProportion (0-1) of sequence that must be covered by high coverage kmers before tagging as redundant"
                    "\n\t\t[--canonical|-c]\t\t\tFlag to ask the program to merge kmers from forward and reverse complement forms (e.g. for DNA-Seq or unstranded RNA-Seq)"
                    "\n\t\t[--filetype|-t (fq|fa; def. fq)]\tWhether the input files are fastq or fasta"
                    "\n\t\t[--outformat|-o (fq|fa; def. fq)]\tWhether you want the output files as fastq or fasta (e.g. for Trinity)"
                    "\n\t\t[--memory_start|-m (integer; def. 1)]\tNumber in Gb of the total memory the program will initially allocate across all threads. The program may request more memory when needed but very small values will cause it to slow down"
                    "\n\t\t[--cpu|-p (int; def 1)]\t\t\tNumber of CPUs that will process the input files, each file is processed sequentially after distributing to the CPUs"
                    "\n\t\t[--verbose|-e]\t\t\t\tEntertain the user"
                    "\n\t\t[--debug|-b]\t\t\t\tAnnoy the developer"
                    "\n\t\t[--version|-v]\t\t\t\tPrint version and exit"
                    "\n\n\n");
}

int parse_arguments(int argc, char *argv[])
{
    cfg.coverage = 0.9;
    cfg.verbose = 0;
    cfg.debug = 0;
    cfg.informat = "fq";
    cfg.is_input_fastq = true;
    cfg.outformat = "fq";
    cfg.is_output_fastq = true;
    cfg.cpus = 1;
    cfg.help = 0;
    cfg.forward_file_count = 0;
    cfg.reverse_file_count = 0;
    cfg.forward_files = NULL;
    cfg.reverse_files = NULL;
    cfg.ksize = 25;
    cfg.depth = 100;
    cfg.depth_per_cpu = cfg.depth / cfg.cpus;
    cfg.memory = 0; // default from INITIAL_CAPACITY - set to be a prime just above 1Gb memory.
    cfg.canonical = 0;
    cfg.version = 0;
    cfg.singleend = 0;

    static struct option long_options[] = {
        {"forward", required_argument, 0, 'f'},
        {"reverse", required_argument, 0, 'r'},
        {"ksize", required_argument, 0, 'k'},
        {"depth", required_argument, 0, 'd'},
        {"coverage", required_argument, 0, 'g'},
        {"filetype", required_argument, 0, 't'},
        {"outformat", required_argument, 0, 'o'},
        {"cpu", required_argument, 0, 'p'},
        {"memory_start", required_argument, 0, 'm'},
        {"debug", required_argument, 0, 'b'},
        {"verbose", no_argument, 0, 'e'},
        {"help", no_argument, 0, 'h'},
        {"canonical", no_argument, 0, 'c'},
        {"version", no_argument, 0, 'v'},
        {"single", no_argument, 0, 's'},
        {"print", no_argument, 0, 'P'},
        {0, 0, 0, 0}};

    int opt;
    int option_index = 0;

    while ((opt = getopt_long(argc, argv, "f:r:k:d:g:t:o:p:m:b:ehcvsP", long_options, &option_index)) != -1)
    {
        switch (opt)
        {
        case 'P':
            cfg.doprint = 1;
            break;
        case 's':
            cfg.singleend = 1;
            break;
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
            print_usage();
            exit(EXIT_SUCCESS);
            break;
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
            cfg.version = 1;
            printf("%d\n", VERSION);
            exit(EXIT_SUCCESS);
            break;
        case 'e':
            cfg.verbose = 1;
            break;
        case 't':
            cfg.informat = optarg;
            if (strcasecmp(cfg.informat, "fa") == 0 || strcasecmp(cfg.informat, "fasta") == 0 || strcasecmp(cfg.informat, "fsa") == 0 || strcasecmp(cfg.informat, "fas") == 0)
            {
                cfg.informat = "fa";
                cfg.is_input_fastq = false;
            }
            else if (strcasecmp(cfg.informat, "fq") == 0 || strcasecmp(cfg.informat, "fastq") == 0 || strcasecmp(cfg.informat, "fsq") == 0)
            {
                cfg.informat = "fq";
                cfg.is_input_fastq = true;
            }
            else
            {
                printf("Input file format must be either fa or fq, not %s\n", cfg.informat);
                return 0;
            }

            break;
        case 'o':
            cfg.outformat = optarg;
            if (strcasecmp(cfg.outformat, "fa") == 0 || strcasecmp(cfg.outformat, "fasta") == 0 || strcasecmp(cfg.outformat, "fsa") == 0 || strcasecmp(cfg.outformat, "fas") == 0)
            {
                cfg.outformat = "fa";
                cfg.is_output_fastq = false;
            }
            else if (strcasecmp(cfg.outformat, "fq") == 0 || strcasecmp(cfg.outformat, "fastq") == 0 || strcasecmp(cfg.outformat, "fsq") == 0)
            {
                cfg.outformat = "fq";
                cfg.is_output_fastq = true;
            }
            else
            {
                printf("Output file format must be either fa or fq, not %s\n", cfg.outformat);
                return 0;
            }
            break;
        default:
            fprintf(stderr, "Unexpected error in option processing\n");
            return 0;
            break;
        }
    }
    if (cfg.verbose)
    {
        printf("\nVERSION: %d, CMD: ", VERSION);
        for (int i = 0; i < argc; i++)
        {
            printf("%s ", argv[i]);
        }
        printf("\n\n");
    }

    // processing after all options parsed:
    cfg.depth_per_cpu = cfg.depth / cfg.cpus;
    // since we're processing each thread chunk separately.
    cfg.initial_hash_size = (cfg.memory > 0) ? memoryGB2capacity(cfg.memory) : INITIAL_CAPACITY;
    float memory_per_cpu = capacity2memory(cfg.initial_hash_size);
    size_t max_slots_for_k = to_power(4, cfg.ksize);
    printf("Initial hash table size set to %'zu (maximum for k=%d is %'zu); memory ~ %'0.2f Gb for each of %d threads (~ %'d Gb total))\n\n", cfg.initial_hash_size, cfg.ksize, max_slots_for_k, memory_per_cpu, cfg.cpus, cfg.memory);

    /////////////////////////////////////////////////////////////////////////
    //////////////      checks       ////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
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
    if (cfg.forward_file_count == 0 || cfg.reverse_file_count == 0)
    {
        fprintf(stderr, "Error: no fwd (%d) or reverse (%d) files provided\n", cfg.forward_file_count, cfg.reverse_file_count);
        return 0;
    }
    if (strcmp(cfg.informat, "fa") == 0 && strcmp(cfg.outformat, "fq") == 0)
    {
        fprintf(stderr, "Error: cannot request an output format of FASTQ when input is FASTA\n");
        return 0;
    }
    if (cfg.singleend == 0 && cfg.forward_file_count != cfg.reverse_file_count)
    {
        fprintf(stderr, "Error: Number of forward (%d) and reverse files (%d) must match\n", cfg.forward_file_count, cfg.reverse_file_count);
        return 0;
    }
    if (cfg.cpus <= 0 || cfg.cpus > MAX_THREADS)
    {
        fprintf(stderr, "Error: CPU count (%d) must be a positive integer and up to %d\n", cfg.cpus, MAX_THREADS);
        return 0;
    }
    if (cfg.ksize < 5 || cfg.ksize > 32)
    {
        fprintf(stderr, "Error: Only kmer sizes (%d) of 5 to 32 are supported\n", cfg.ksize);
        return 0;
    }
    if (cfg.coverage > 1 || cfg.coverage < 0.001)
    {
        fprintf(stderr, "Error: Coverage (%3.f) is the proportion of the sequence covered by high kmers and must be between 0 and 1\n", cfg.coverage);
        return 0;
    }
    if (cfg.depth < 2)
    {
        fprintf(stderr, "Error: Depth (%d) is the number of times a kmer needs to be found before being flagged as high coverage, it must be above 1\n", cfg.depth);
        return 0;
    }
    if (cfg.depth_per_cpu < 2)
    {
        fprintf(stderr, "Error: Depth (%d) must be at least 2 x number of CPUs (for performance reasons; but this version of the program is written to normalise to 50+\n", cfg.depth);
        return 0;
    }

    if (cfg.initial_hash_size < 100000)
    {
        fprintf(stderr, "Error: initial kmer table size is too small (%'zu), it should be set to at least 100000 (or leave empty for default)\n", cfg.initial_hash_size);
        return 0;
    }

    return 1;
}

void process_forward_files(char *first_file, char **additional_files)
{
    // the first file
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

    // additional files
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
    // i need to fix this one day
    cfg.reverse_files = realloc(cfg.reverse_files, cfg.forward_file_count * sizeof(char *));

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

char *create_output_filename(const char *basename, int k, int norm_depth, int thread, char *suffix)
{
    size_t len = strlen(basename) + 50;
    char *output_filename = malloc(len);
    if (!output_filename)
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    if (thread >= 0)
        snprintf(output_filename, len, "%s.k%d_norm%d_thread%d.%s", basename, k, norm_depth, thread, suffix);
    else
        snprintf(output_filename, len, "%s.k%d_norm%d.%s", basename, k, norm_depth, suffix);

    return output_filename;
}

void fastq_to_fasta(char (*fastq_record)[4][MAX_LINE_LENGTH], char *fasta_output, bool is_forward)
{
    const char *suffix = (is_forward) ? "/1" : "/2";

    char header[MAX_LINE_LENGTH];
    char seq[MAX_LINE_LENGTH];

    fasta_output[0] = '\0';

    strcpy(header, (*fastq_record)[0]);
    strcpy(seq, (*fastq_record)[1]);

    header[0] = '>';

    size_t len = strlen(header);
    if (len < 2 || strcmp(header + len - 2, suffix) != 0)
    {
        strcat(header, suffix);
    }

    strcat(fasta_output, header);
    strcat(fasta_output, "\n");
    strcat(fasta_output, seq);
    strcat(fasta_output, "\n");
}

////////////////////////////////////////////////////////////////
////////////////HASH TABLE FUNCTIONS ///////////////////////////
////////////////////////////////////////////////////////////////
void destroy_hash_table(hash_table_t *ht)
{
    if (ht && ht->kmer)
    {
        free(ht->kmer);
        ht->kmer = NULL;
    }
}

void init_hash_table(hash_table_t *ht)
{
    if (ht && ht->kmer && ht->used > 0)
    {
        printf("Existing hash table with %zu entries found (shouldn't happen)\n", ht->used);
        exit(EXIT_FAILURE);
    }

    ht->kmer = calloc(cfg.initial_hash_size, sizeof(kmer_t));
    if (!ht->kmer)
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
    copy->kmer = calloc(copy->capacity, sizeof(kmer_t));
    if (!copy->kmer)
    {
        fprintf(stderr, "Memory allocation failed\n");
        free(copy);
        exit(EXIT_FAILURE);
    }
    memcpy(copy->kmer, source->kmer, copy->capacity * sizeof(kmer_t));
    return copy;
}

size_t store_kmer(hash_table_t *hash_table, uint64_t hash, int thread_id, bool do_init, char *kmer)
{
    // this happens locally on a thread specific hash_table so no requirement to do any locking

    if (hash_table->used >= hash_table->capacity * TABLE_LOAD_FACTOR)
        expand_local_hash_table(hash_table, 0, thread_id);

    size_t index = hash % hash_table->capacity;

    if (index < 0)
    {
        fprintf(stderr, "ERROR1: This shouldnt have happened, index %zu hash %zu capacity %zu\n", index, hash, hash_table->capacity);
        exit(EXIT_FAILURE);
    }

    if (cfg.debug > 2)
        printf("DEBUG: Kmer hash: %lu, Count: %d\n", hash, hash_table->kmer[index].count);

    // if entry is null and available:
    if (hash_table->kmer[index].hash == 0)
    {
        if (cfg.debug > 3)
        {
            char kmer_str[cfg.ksize + 1];
            decode_kmer_plain(hash, cfg.ksize, kmer_str);
            printf("Thread %d: new kmer %s derived from %s, hash %zu (existing: %zu) at index %zu. Existing count is %d, and hash capacity is %zu and used size %zu", thread_id, kmer_str, kmer, hash, hash_table->kmer[index].hash, index, hash_table->kmer[index].count, hash_table->capacity, hash_table->used);

            if (strcmp(kmer_str, kmer) != 0)
            {
                printf("\n Thes two kmers are not the same!\n");
                exit(EXIT_FAILURE);
            }
        }
        hash_table->kmer[index].hash = hash;
        hash_table->kmer[index].count = (do_init == true) ? (int)0 : 1;
        // new kmer inserted, so increase size used.
        hash_table->used++;

        if (cfg.debug > 3)
            printf(" new count is %d, and new used size %zu\n", hash_table->kmer[index].count, hash_table->used);

        return index;
    }
    else if (hash_table->kmer[index].hash == hash)
    {
        if (cfg.debug > 3)
        {
            char kmer_str[cfg.ksize + 1];
            decode_kmer_plain(hash, cfg.ksize, kmer_str);
            printf("Thread %d: existing kmer %s derived from %s, hash %'zu (existing: %'zu) at index %zu. Existing count is %d, and hash capacity is %'zu and used size %'zu",
                   thread_id, kmer_str, kmer, hash, hash_table->kmer[index].hash, index, hash_table->kmer[index].count, hash_table->capacity, hash_table->used);

            // if (thread_id < 0 && hash_table->kmer[index].count != 0)
            // {
            //     printf(" This should have happened\n");
            //     exit(EXIT_FAILURE);
            // }

            if (strcmp(kmer_str, kmer) != 0)
            {
                printf("\n These two kmers are not the same!\n");
                exit(EXIT_FAILURE);
            }
        }

        if (do_init == true)
        {
            if (cfg.debug > 3)
            {
                printf(" seed initialisation so not updating count.\n");
            }
            return index;
        }
        else
            hash_table->kmer[index].count++;
        if (cfg.debug > 3)
        {
            printf(" new count is %d, and new used size %'zu\n", hash_table->kmer[index].count, hash_table->used);
        }
    }
    else // collision
    {

        size_t original_index = index;
        int collisions = 0;

        while (hash_table->kmer[index].hash != 0 && hash_table->kmer[index].hash != hash)
        {
            collisions++;
            if (cfg.debug > 4)
                printf("Thread %d: hash %'zu collision consecutive number %d, index: %'zu -> %'zu, capacity %'zu\n", thread_id, hash, collisions, original_index, index, hash_table->capacity);

            if ((collisions / hash_table->capacity) > (hash_table->capacity * 0.5))
            {
                if (cfg.debug)
                    printf("Thread %d: Collisions more than 10%% of table capacity, expanding table...\n", thread_id);
                expand_local_hash_table(hash_table, 0, thread_id);
            }

            index = (index + 1) % hash_table->capacity;
            if (index < 0)
            {
                fprintf(stderr, "ERROR1: This shouldnt have happened, index %zu hash %zu capacity %zu\n", index, hash, hash_table->capacity);
                exit(EXIT_FAILURE);
            }

            if (do_init == false)
                hash_table->kmer[index].count++;
            else
                hash_table->kmer[index].count = (int)0;
        }
    }

    if (cfg.debug > 2)
        printf("DEBUG: New Kmer hash: %lu, Count: %d\n", hash, hash_table->kmer[index].count);
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
        if (ht->kmer[i].hash != 0)
        {
            // kmer is stored at an index based on its capacity so when that increases, indexes change
            size_t new_index = ht->kmer[i].hash % new_capacity;
            while (new_entries[new_index].hash != 0)
            {
                new_index = (new_index + 1) % new_capacity;
            }
            new_entries[new_index] = ht->kmer[i];
            new_hash_size++;
        }
    }

    if (ht->kmer)
        free(ht->kmer);

    ht->kmer = new_entries;
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
// is this good enough? stop char is > or @ and @ can occur in quality seqs.
// char *find_next_record_start(char *ptr, char *end_ptr, char stop_char)
// {
//     while (ptr < end_ptr && *ptr != stop_char)
//     {
//         ptr = strchr(ptr, '\n');
//         if (ptr)
//             ptr++;
//         else
//             break;
//     }
//     return (ptr < end_ptr) ? ptr : NULL;
// }

size_t find_thread_exact_end(char *data, size_t start_pos, size_t end_pos)
{

    // fasta
    if (cfg.is_input_fastq == false)
    {
        for (size_t i = end_pos; i > start_pos; i--) // loop backwards
            if (data[i] == '>')
                return i - 1;
    }
    else
    {
        // fastq
        int line_break_count = 0;
        bool header_qual_found = false;
        for (size_t i = end_pos; i > start_pos; i--) // loop backwards
        {
            // printf("1:checking start position at %zu, line_break_count %d, line is %c\n", i, line_break_count, data[i]);
            if (data[i] == '\n') // when we have a new line, check the next character.
            {
                line_break_count++;
                if (data[i + 1] == '+')
                    header_qual_found = true;
                else if (header_qual_found && data[i + 1] == '@')
                    return i;
                if (line_break_count == 7)
                {
                    printf("ERROR: after 7 lines, I couldn't find the + and @ headers near this chunk %zu\n", i);
                    exit(EXIT_FAILURE);
                }
            }
        }
    }

    printf("ERROR: i couldn't find the start of sequence before this chunk end %zu\n", end_pos);
    exit(EXIT_FAILURE);
    return -1;
}

// this has a bug, we don't account for forward and reverse files having different sizes.
// we should just do line counting...??!?!?
void calculate_thread_positions(char *data, size_t total_file_size, int thread_count, size_t **thread_starts, size_t **thread_ends)
{
    size_t approx_chunk_size = total_file_size / thread_count;
    size_t approximate_end = approx_chunk_size - (MAX_LINE_LENGTH * 4);

    //    printf("file size is %'zu thread count is %d\n", total_file_size, thread_count);

    (*thread_starts)[0] = 0;
    (*thread_ends)[0] = find_thread_exact_end(data, 0, approximate_end);
    (*thread_ends)[thread_count - 1] = total_file_size - 1; // last thread, end of file

    for (int thread_number = 1; thread_number < thread_count; thread_number++)
    {
        size_t start_pos = (*thread_ends)[thread_number - 1] + 1; // start is end of previous one
        size_t end_pos = start_pos + approximate_end;
        // printf("thread %d search space: start_pos %'zu end_pos %'zu\n", i, start_pos, end_pos);

        (*thread_ends)[thread_number] = find_thread_exact_end(data, start_pos, end_pos);
        // printf("thread %d start set to %'zu end to %'zu\n", i, (*thread_starts)[i],(*thread_ends)[i]);
        if (thread_number < thread_count - 1)
            (*thread_starts)[thread_number + 1] = (*thread_ends)[thread_number] + 1;
    }
}

// a different and better approach would be to submit a thread for X records rather than parse the file and split across. another day
void calculate_thread_positions_from_records(char *data, size_t total_file_size, int thread_count, size_t **thread_starts, size_t **thread_ends, size_t records)
{

    size_t records_per_thread = records / thread_count; // last thread will have a bit more

    if (thread_count < 2 || records_per_thread < 1 || total_file_size < 1)
        return;

    int max_lines = cfg.is_input_fastq == true ? records_per_thread * 4 : records_per_thread * 2;

    // thread 0 is 0
    (*thread_starts)[0] = 0;
    (*thread_ends)[thread_count - 1] = total_file_size - 1; // last thread, end of file

    // don't need last thread
    for (int thread_number = 0; thread_number < thread_count - 1; thread_number++)
    {
        size_t line_count = 0;
        for (size_t i = (*thread_starts)[thread_number]; i < total_file_size; i++)
        {
            if (data[i] == '\n')
            {
                line_count++;
                if (line_count == max_lines)
                {
                    (*thread_ends)[thread_number] = i;

                    if (thread_number < thread_count - 1)
                        (*thread_starts)[thread_number + 1] = i + 1;

                    break;
                }
            }
        }
    }
}

size_t count_records_seqfile(char *data, size_t size)
{
    size_t line_count = 0;
    for (size_t i = 0; i < size; ++i)
    {
        if (data[i] == '\n')
        {
            line_count++;
        }
    }
    // Add 1 for the last line if it doesn't end with a newline
    if (size > 0 && data[size - 1] != '\n')
        line_count++;

    if (cfg.is_input_fastq == true)
        return line_count / 4;
    else
        return line_count / 2;
}

void seed_kmer_hash(char *file, int records_to_seed, hash_table_t *ht, int K, int depth_per_cpu)
{
    if (cfg.verbose)
        printf("Seeding hash table with up to %'d records from file %s\n", records_to_seed, file);

    mmap_file_t mmapped = mmap_file(file);

    int lines_to_read = cfg.is_input_fastq == true ? 4 : 2;
    int record_count = 0, line_count = 0;
    char record[lines_to_read][MAX_LINE_LENGTH];
    size_t char_index = 0;

    for (size_t character = 0; character < mmapped.size; ++character)
    {
        record[line_count][char_index++] = mmapped.data[character];

        if (mmapped.data[character] == '\n')
        {
            record[line_count][char_index - 1] = '\0'; // replace newline with \0
            char_index = 0;
            line_count++;

            if (line_count == lines_to_read)
            {

                if (strlen(record[1]) > K)
                {
                    if (is_valid_sequence_single(record[1], K) == false)
                        exit(EXIT_FAILURE);

                    sequence_to_hash_zero(record[1], ht, K, depth_per_cpu);
                    record_count++;
                    if (record_count == records_to_seed)
                    {
                        if (cfg.verbose)
                            printf("%'d records with %'zu kmers seeded\n", record_count, ht->used);
                        munmap_file(&mmapped);
                        return;
                    }
                }
                line_count = 0;
                for (int i = 0; i < lines_to_read; ++i)
                    record[i][0] = '\0';
            }
        }
    }
    munmap_file(&mmapped);
    if (cfg.verbose)
        printf("%'d records with %'zu kmers seeded\n", record_count, ht->used);

    return;
}

// size_t get_position_from_sequence_record(char *data, size_t size, size_t record_number)
// {
//     size_t line_number = 0;

//     if (cfg.is_input_fastq == true)
//         line_number = record_number * 4;
//     else
//         line_number = record_number * 2;

//     size_t current_line = 1;
//     size_t current_position = 0;

//     for (size_t i = 0; i < size; ++i)
//     {
//         if (data[i] == '\n')
//         {
//             current_line++;
//             current_position = i;
//             if (current_line == line_number)
//             {
//                 printf("Line %'zu found at position %zu, data: %c", line_number, current_position, data[i]);
//                 return current_position + 1;
//             }
//         }
//     }
//     printf("ERROR: Line number %zu not found in file\n", line_number);
//     exit(EXIT_FAILURE);
//     return (size_t)-1;
// }
bool is_valid_sequence_single(const char *fwd_seq, int K)
{
    replacestr(fwd_seq, "N", "A");
    int fwd_seq_len = strlen(fwd_seq);
    if (fwd_seq_len < K)
    {
        // if (cfg.debug > 1)
        //     printf("Thread %d - Sequence FWD skipped due to length shorter than K\n%s\n\n", thread_id, fwd_seq);

        return false;
    }

    if (valid_dna(fwd_seq) == false)
    {
        fprintf(stderr, "FATAL: FWD sequence does not appear to be a DNA sequence\n%s\n\n", fwd_seq);
        exit(EXIT_FAILURE);
    }

    return true;
}
bool is_valid_sequence_pair(const char *fwd_seq, const char *rev_seq, int K)
{
    replacestr(fwd_seq, "N", "A");
    replacestr(rev_seq, "N", "A");
    int fwd_seq_len = strlen(fwd_seq);
    int rev_seq_len = strlen(rev_seq);
    if (fwd_seq_len < K)
    {
        // if (cfg.debug > 1)
        //     printf("Thread %d - Sequence FWD skipped due to length shorter than K\n%s\n\n", thread_id, fwd_seq);

        return false;
    }
    if (rev_seq_len < K)
    {
        // if (cfg.debug > 1)
        //     printf("Thread %d - Sequence REV skipped due to length shorter than K\n%s\n\n", thread_id, rev_seq);

        return false;
    }

    if (valid_dna(fwd_seq) == false)
    {
        fprintf(stderr, "FATAL: FWD sequence does not appear to be a DNA sequence\n%s\n\n", fwd_seq);
        exit(EXIT_FAILURE);
    }
    if (valid_dna(rev_seq) == false)
    {
        fprintf(stderr, "FATAL: REV sequence does not appear to be a DNA sequence\n%s\n\n", rev_seq);
        exit(EXIT_FAILURE);
    }

    return true;
}

void sequence_to_hash(const char *seq, hash_table_t *hash_table, int *seq_high_count_kmers, int *total_seq_kmers, int K, int depth_per_cpu, int thread_id)
{
    // if (cfg.debug > 3)
    // printf("Kmerising Sequence %s\n", seq);

    for (int i = 0; i <= strlen(seq) - K; i++)
    {
        char kmer[K + 1];
        strncpy(kmer, seq + i, K);
        kmer[K] = '\0';

        uint64_t hash = 0;

        if (cfg.canonical == 1)
        {
            const char *canonical_kmer = get_canonical_kmer(kmer, K);
            hash = encode_kmer_plain(canonical_kmer, K);
        }
        else
        {
            hash = encode_kmer_plain(kmer, K);
        }

        // this actually means hashing failed so it shouldn't occur but check so we don't get segfault
        if (hash == 0)
            continue;

        (*total_seq_kmers)++;
        size_t index = store_kmer(hash_table, hash, thread_id, false, kmer);

        if (index < 0)
        {
            fprintf(stderr, "ERROR3: This shouldnt have happened, index %zu hash %zu capacity %zu\n", index, hash, hash_table->capacity);
            exit(EXIT_FAILURE);
        }
        if (hash_table->kmer[index].count >= depth_per_cpu)
        {
            (*seq_high_count_kmers)++;
        }
    }
}

void sequence_to_hash_zero(const char *seq, hash_table_t *hash_table, int K, int depth_per_cpu)
{
    // if (cfg.debug > 3)
    // printf("Kmerising Sequence %s\n", seq);

    for (int i = 0; i <= strlen(seq) - K; i++)
    {
        char kmer[K + 1];
        strncpy(kmer, seq + i, K);
        kmer[K] = '\0';
        // printf("length is %zu (%s)\n", strlen(kmer), kmer);

        uint64_t hash = 0;

        if (cfg.canonical == 1)
        {
            const char *canonical_kmer = get_canonical_kmer(kmer, K);
            hash = encode_kmer_plain(canonical_kmer, K);
        }
        else
        {
            hash = encode_kmer_plain(kmer, K);
        }

        // this actually means hashing failed so it shouldn't occur but check so we don't get segfault
        if (hash == 0)
            continue;

        size_t index = store_kmer(hash_table, hash, -1, true, kmer);

        if (index < 0)
        {
            fprintf(stderr, "ERROR3: This shouldnt have happened, index %zu hash %zu capacity %zu\n", index, hash, hash_table->capacity);
            exit(EXIT_FAILURE);
        }
    }
}

bool process_sequence_single(const char *fwd_seq, hash_table_t *hash_table, int *fwd_seq_high_count_kmers,
                             int *fwd_total_seq_kmers, int thread_id, int K, int depth_per_cpu)
{
    if (is_valid_sequence_single(fwd_seq, K) == false)
        return false;

    *fwd_seq_high_count_kmers = 0;
    *fwd_total_seq_kmers = 0;
    sequence_to_hash(fwd_seq, hash_table, fwd_seq_high_count_kmers, fwd_total_seq_kmers, K, depth_per_cpu, thread_id);

    return true;
}
bool process_sequence_pair(const char *fwd_seq, const char *rev_seq, hash_table_t *hash_table, int *fwd_seq_high_count_kmers,
                           int *fwd_total_seq_kmers, int *rev_seq_high_count_kmers, int *rev_total_seq_kmers, int thread_id, int K, int depth_per_cpu)
{
    if (is_valid_sequence_pair(fwd_seq, rev_seq, K) == false)
        return false;

    *fwd_seq_high_count_kmers = 0;
    *fwd_total_seq_kmers = 0;
    sequence_to_hash(fwd_seq, hash_table, fwd_seq_high_count_kmers, fwd_total_seq_kmers, K, depth_per_cpu, thread_id);

    *rev_seq_high_count_kmers = 0;
    *rev_total_seq_kmers = 0;
    sequence_to_hash(rev_seq, hash_table, rev_seq_high_count_kmers, rev_total_seq_kmers, K, depth_per_cpu, thread_id);

    return true;
}

void *process_thread_chunk_paired(void *arg)
{
    thread_data_t *data = (thread_data_t *)arg;
    // generic
    int lines_to_read = cfg.is_input_fastq == true ? 4 : 2;
    hash_table_t *local_hash_table = data->hash_table;
    int thread_id = data->thread_id;
    FILE *output_forward = data->thread_output_forward;
    FILE *output_reverse = data->thread_output_reverse;
    // so we don't keep accessing it.
    int K = cfg.ksize;
    int depth_per_cpu = cfg.depth_per_cpu;
    int debug = cfg.debug;
    int verbose = cfg.verbose;
    int singleend = cfg.singleend;
    bool is_input_fastq = cfg.is_input_fastq;
    bool is_output_fastq = cfg.is_output_fastq;
    float coverage = cfg.coverage;

    // input file
    char *forward_pointer = data->forward_data + data->forward_file_start;
    char *reverse_pointer = data->reverse_data + data->reverse_file_start;

    // thread stats
    size_t total_kmers = data->hash_table->used;
    int processed_count = 0, printed_count = 0, skipped_count = 0, prev_printed_count = data->printed_count, prev_skipped_count = data->skipped_count;
    size_t prev_total_kmers = total_kmers, sequences_last_processed = data->processed_count;
    double prev_rate = 0;
    size_t sequences_processed = 0;
    time_t current_time = time(NULL);
    data->last_report_time = current_time;
    data->last_report_count = data->processed_count;

    if (cfg.debug > 1)
        printf("Thread %d started; processing paired files, %d lines per record\n", thread_id, lines_to_read);

    // all checks in here may slow down things
    while (forward_pointer < data->forward_data + data->forward_file_end &&
           reverse_pointer < data->reverse_data + data->reverse_file_end)
    {

        bool valid_record = true;
        int seq_high_count_kmers_forward = 0, total_seq_kmers_forward = 0;
        int seq_high_count_kmers_reverse = 0, total_seq_kmers_reverse = 0;

        char forward_record[lines_to_read][MAX_LINE_LENGTH], reverse_record[lines_to_read][MAX_LINE_LENGTH];

        // read pair
        for (int i = 0; i < lines_to_read; i++)
        {
            // store data into record, update pointer
            forward_pointer = read_line(forward_pointer, forward_record[i], MAX_LINE_LENGTH);
            reverse_pointer = read_line(reverse_pointer, reverse_record[i], MAX_LINE_LENGTH);

            if (!forward_pointer || !reverse_pointer)
            {
                valid_record = false;
                break;
            }
        }
        if (!valid_record)
            break;

        bool process_success = process_sequence_pair(forward_record[1], reverse_record[1], local_hash_table, &seq_high_count_kmers_forward, &total_seq_kmers_forward, &seq_high_count_kmers_reverse, &total_seq_kmers_reverse, thread_id, K, depth_per_cpu);
        if (process_success == false)
            continue;

        data->processed_count++;

        // note, we process and filter pairs; it's the average redundancy in the pair not each sequence.
        // i think this is appropriate but it could be replaced with an AND clause (but not an OR clause)

        // int seq_high_count_kmers_pair = seq_high_count_kmers_forward + seq_high_count_kmers_reverse;
        // int total_seq_kmers_pair = total_seq_kmers_forward + total_seq_kmers_reverse;
        // float high_count_ratio_pair = total_seq_kmers_pair > 0 ? (float)seq_high_count_kmers_pair / total_seq_kmers_pair : 0;
        float high_count_ratio_fwd = total_seq_kmers_forward > 0 ? (float)seq_high_count_kmers_forward / total_seq_kmers_forward : 0;
        float high_count_ratio_rev = total_seq_kmers_reverse > 0 ? (float)seq_high_count_kmers_reverse / total_seq_kmers_reverse : 0;

        bool is_printed = false;
        // options
        if (high_count_ratio_fwd <= coverage && high_count_ratio_rev <= coverage)
        // if (high_count_ratio_pair <= coverage) // this skips fewer sequences, which is a tad surprising.
        {
            if (is_input_fastq == true && is_output_fastq == false)
            {
                char forward_record_fasta[MAX_LINE_LENGTH * 2] = {0};
                fastq_to_fasta(&forward_record, forward_record_fasta, true);

                char reverse_record_fasta[MAX_LINE_LENGTH * 2] = {0};
                fastq_to_fasta(&reverse_record, reverse_record_fasta, false);
                fprintf(output_forward, "%s", forward_record_fasta);
                fprintf(output_reverse, "%s", reverse_record_fasta);
            }
            else
            {
                for (int i = 0; i < lines_to_read; i++)
                {
                    fprintf(output_forward, "%s\n", forward_record[i]);
                    fprintf(output_reverse, "%s\n", reverse_record[i]);
                }
            }
            data->printed_count++;
            is_printed = true;
        }
        else
        {
            data->skipped_count++;
            is_printed = false;
        }

        // debug report (this check isn't a massive hit)
        if (debug > 1)
        {

            if (is_printed == true)
            {
                printf("Thread %d - Sequence pair %'zu PRINTED: ", thread_id, data->processed_count);
            }
            else
            {
                printf("Thread %d - Sequence pair %'zu SKIPPED: ", thread_id, data->processed_count);
            }
            printf("High (%d) count kmers: F:%d;R:%d, Total unique kmers: F:%d;R:%d High count ratio: F:%.2f;R:%.2f\n",
                   depth_per_cpu,
                   seq_high_count_kmers_forward, seq_high_count_kmers_reverse,
                   total_seq_kmers_forward, total_seq_kmers_reverse,
                   high_count_ratio_fwd, high_count_ratio_rev);
        }
        sequences_processed++;

        current_time = time(NULL);
        if (difftime(current_time, data->last_report_time) >= REPORTING_INTERVAL)
        {
            if (debug > 1)
                printf("reporting after %d seconds\n", REPORTING_INTERVAL);

            double elapsed_time = difftime(current_time, data->last_report_time);
            sequences_last_processed = data->processed_count - data->last_report_count;
            double rate = sequences_last_processed / elapsed_time;
            total_kmers = local_hash_table->used;

            float printed_improvement = (prev_printed_count == 0) ? 0 : (float)(data->printed_count - prev_printed_count) / prev_printed_count;
            float skipped_improvement = (prev_skipped_count == 0) ? 0 : (float)(data->skipped_count - prev_skipped_count) / prev_skipped_count;
            float prev_rate_improvement = (prev_rate == 0) ? 0 : (float)(rate - prev_rate) / prev_rate;
            float kmer_improvement = (prev_total_kmers == 0) ? 0 : (float)(total_kmers - prev_total_kmers) / prev_total_kmers;

            if (verbose)
                printf("Thread %d - Processing rate: %'.0f (%+.2f%%) sequences/s, processed %'zu pairs, printed: %'zu (%+.2f%%), skipped: %'zu (%+.2f%%), Unique kmers (all sequences; this thread): %'zu (%+.2f%%)\n",
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
    }

    // chunk finished processing

    if (verbose)
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

    printf("Thread %d - Processing rate: %'.0f (%+.2f%%) sequences/s, processed %'zu pairs, printed: %'zu (%+.2f%%), skipped: %'zu (%+.2f%%), Unique kmers (all sequences; this thread): %'zu (%+.2f%%)\n",
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

int multithreaded_process_files_paired(thread_data_t *thread_data, mmap_file_t *forward_mmap, mmap_file_t *reverse_mmap)
{
    char stop_char = cfg.is_input_fastq == false ? '>' : '@';

    size_t max_total_kmers_in_threads = 0;
    pthread_t threads[cfg.cpus];

    // let's get start and ends for the mmapped file.
    size_t *forward_thread_starts = NULL;
    size_t *forward_thread_ends = NULL;
    size_t *reverse_thread_starts = NULL;
    size_t *reverse_thread_ends = NULL;
    forward_thread_starts = calloc(cfg.cpus, sizeof(size_t));
    forward_thread_ends = calloc(cfg.cpus, sizeof(size_t));
    reverse_thread_starts = calloc(cfg.cpus, sizeof(size_t));
    reverse_thread_ends = calloc(cfg.cpus, sizeof(size_t));

    if (!forward_thread_starts || !forward_thread_ends || !reverse_thread_starts || !reverse_thread_ends)
    {
        perror("Error allocating memory for thread positions");
        exit(EXIT_FAILURE);
    }

    // todo: bug, nb this will not work because forward and reverse can have different file sizes.
    if (cfg.cpus == 1)
    {
        if (cfg.verbose)
            printf("Single thread mode\n");

        forward_thread_ends[0] = forward_mmap->size - 1;
        reverse_thread_ends[0] = reverse_mmap->size - 1;
    }
    else
    {

        if (forward_mmap->size == reverse_mmap->size)
        {
            if (cfg.verbose)
                printf("The forward and reverse files have the same file size, assuming same number of records!\n");

            calculate_thread_positions(forward_mmap->data, forward_mmap->size, cfg.cpus, &forward_thread_starts, &forward_thread_ends);
            calculate_thread_positions(reverse_mmap->data, reverse_mmap->size, cfg.cpus, &reverse_thread_starts, &reverse_thread_ends);
        }
        else
        {
            if (cfg.verbose)
                printf("The forward (%'zu) and reverse files (%'zu) have different file size, so calculating split amongst threads is slower, hold on...\n",
                       forward_mmap->size, reverse_mmap->size);

            size_t total_records = 0;
            total_records = count_records_seqfile(forward_mmap->data, forward_mmap->size);

            if (cfg.debug > 0)
                printf("forward file has %'zu records\n", total_records);

            calculate_thread_positions_from_records(forward_mmap->data, forward_mmap->size, cfg.cpus, &forward_thread_starts, &forward_thread_ends, total_records);
            calculate_thread_positions_from_records(reverse_mmap->data, reverse_mmap->size, cfg.cpus, &reverse_thread_starts, &reverse_thread_ends, total_records);

            for (int thread_number = 0; thread_number < cfg.cpus; thread_number++)
            {
                if (cfg.debug > 0)
                {
                    printf("Thread %d: fwd record starts at %'zu (%c) and ends at %'zu (newline after %c)\n", thread_number, forward_thread_starts[thread_number], forward_mmap->data[forward_thread_starts[thread_number]], forward_thread_ends[thread_number], forward_mmap->data[forward_thread_ends[thread_number] - 1]);
                    printf("Thread %d: rev record starts at %'zu (%c) and ends at %'zu (newline after %c)\n", thread_number, reverse_thread_starts[thread_number], reverse_mmap->data[reverse_thread_starts[thread_number]], reverse_thread_ends[thread_number], reverse_mmap->data[reverse_thread_ends[thread_number] - 1]);
                }
            }
        }
    }

    for (int thread_number = 0; thread_number < cfg.cpus; thread_number++)
    {

        // store memory mapped file starts and stops, and data structure.
        thread_data[thread_number].forward_file_start = forward_thread_starts[thread_number];
        thread_data[thread_number].forward_file_end = forward_thread_ends[thread_number];
        thread_data[thread_number].reverse_file_start = reverse_thread_starts[thread_number];
        thread_data[thread_number].reverse_file_end = reverse_thread_ends[thread_number];
        thread_data[thread_number].forward_data = forward_mmap->data;
        thread_data[thread_number].reverse_data = reverse_mmap->data;

        if (!thread_data[thread_number].thread_output_forward || !thread_data[thread_number].thread_output_reverse)
        {
            fprintf(stderr, "Error opening thread-specific output files for thread %d, %s and %s\n", thread_number, thread_data[thread_number].thread_output_forward_filename, thread_data[thread_number].thread_output_reverse_filename);
            exit(EXIT_FAILURE);
        }

        // printf("forward/reverse size %zu/%zu chunk %s/%s\n", thread_data[i].forward_size, thread_data[i].reverse_size, forward_chunk_end, reverse_chunk_end);

        if (cfg.debug > 0)
            printf("Starting thread %d\n", thread_number);

        if (pthread_create(&threads[thread_number], NULL, process_thread_chunk_paired, &thread_data[thread_number]) != 0)
        {

            fprintf(stderr, "Error creating thread %d\n", thread_number);
            for (int j = 0; j < thread_number; j++)
            {
                pthread_cancel(threads[j]);
                pthread_join(threads[j], NULL);
            }
            free(forward_thread_starts);
            free(forward_thread_ends);
            free(reverse_thread_starts);
            free(reverse_thread_ends);
            return 1;
        }

        sleep(1);
    }

    // close
    for (int i = 0; i < cfg.cpus; i++)
    {
        if (pthread_join(threads[i], NULL) != 0)
        {
            fprintf(stderr, "Error joining thread %d\n", i);
            free(forward_thread_starts);
            free(forward_thread_ends);
            free(reverse_thread_starts);
            free(reverse_thread_ends);
            return 1;
        }
    }

    // final report
    size_t total_processed = 0, total_printed = 0, total_skipped = 0;
    for (int thread_number = 0; thread_number < cfg.cpus; thread_number++)
    {
        total_processed += thread_data[thread_number].processed_count;
        total_printed += thread_data[thread_number].printed_count;
        total_skipped += thread_data[thread_number].skipped_count;
        max_total_kmers_in_threads = (max_total_kmers_in_threads < thread_data[thread_number].hash_table->used) ? thread_data[thread_number].hash_table->used : max_total_kmers_in_threads;

        if (cfg.doprint)
            print_kmer_table(thread_data[thread_number].hash_table, "", true, thread_number);
    }
    reporting.total_processed = total_processed;
    reporting.total_printed = total_printed;
    reporting.total_skipped = total_skipped;
    reporting.files_processed++;
    reporting.max_total_kmers = (reporting.max_total_kmers < max_total_kmers_in_threads) ? max_total_kmers_in_threads : reporting.max_total_kmers;

    printf("Cumulative file statistics: Processed %'zu, Printed %'zu, Skipped %'zu, Cumulative Max Unique Kmers in a thread: %'zu\n",
           total_processed, total_printed, total_skipped, max_total_kmers_in_threads);

    free(forward_thread_starts);
    free(forward_thread_ends);
    free(reverse_thread_starts);
    free(reverse_thread_ends);

    return 0;
}
void *process_thread_chunk_single(void *arg)
{
    thread_data_t *data = (thread_data_t *)arg;
    // generic
    int lines_to_read = cfg.is_input_fastq == true ? 4 : 2;
    hash_table_t *local_hash_table = data->hash_table;
    int thread_id = data->thread_id;
    FILE *output_forward = data->thread_output_forward;
    // so we don't keep accessing it.
    int K = cfg.ksize;
    int depth_per_cpu = cfg.depth_per_cpu;
    int debug = cfg.debug;
    int verbose = cfg.verbose;
    bool is_input_fastq = cfg.is_input_fastq;
    bool is_output_fastq = cfg.is_output_fastq;
    float coverage = cfg.coverage;

    // input file
    char *forward_pointer = data->forward_data + data->forward_file_start;

    // thread stats
    size_t total_kmers = data->hash_table->used;
    int processed_count = 0, printed_count = 0, skipped_count = 0, prev_printed_count = data->printed_count, prev_skipped_count = data->skipped_count;
    size_t prev_total_kmers = total_kmers, sequences_last_processed = data->processed_count;
    double prev_rate = 0;
    size_t sequences_processed = 0;
    time_t current_time = time(NULL);
    data->last_report_time = current_time;
    data->last_report_count = data->processed_count;

    if (cfg.debug > 1)
        printf("Thread %d started; processing single end file, %d lines per record\n", thread_id, lines_to_read);

    // all checks in here may slow down things
    while (forward_pointer < data->forward_data + data->forward_file_end)
    {

        bool valid_record = true;
        int seq_high_count_kmers_forward = 0, total_seq_kmers_forward = 0;

        char forward_record[lines_to_read][MAX_LINE_LENGTH], reverse_record[lines_to_read][MAX_LINE_LENGTH];

        // read pair
        for (int i = 0; i < lines_to_read; i++)
        {
            // store data into record, update pointer
            forward_pointer = read_line(forward_pointer, forward_record[i], MAX_LINE_LENGTH);

            if (!forward_pointer)
            {
                valid_record = false;
                break;
            }
        }
        if (!valid_record)
            break;
        bool process_success = process_sequence_single(forward_record[1], local_hash_table, &seq_high_count_kmers_forward, &total_seq_kmers_forward, thread_id, K, depth_per_cpu);
        if (process_success == false)
            continue;

        data->processed_count++;

        // note, we process and filter pairs; it's the average redundancy in the pair not each sequence.
        // i think this is appropriate but it could be replaced with an AND clause (but not an OR clause)

        // int seq_high_count_kmers_pair = seq_high_count_kmers_forward + seq_high_count_kmers_reverse;
        // int total_seq_kmers_pair = total_seq_kmers_forward + total_seq_kmers_reverse;
        // float high_count_ratio_pair = total_seq_kmers_pair > 0 ? (float)seq_high_count_kmers_pair / total_seq_kmers_pair : 0;
        float high_count_ratio_fwd = total_seq_kmers_forward > 0 ? (float)seq_high_count_kmers_forward / total_seq_kmers_forward : 0;

        bool is_printed = false;
        // options
        if (high_count_ratio_fwd <= coverage)
        // if (high_count_ratio_pair <= coverage) // this skips fewer sequences, which is a tad surprising.
        {
            if (is_input_fastq == true && is_output_fastq == false)
            {
                char forward_record_fasta[MAX_LINE_LENGTH * 2] = {0};
                fastq_to_fasta(&forward_record, forward_record_fasta, true);

                char reverse_record_fasta[MAX_LINE_LENGTH * 2] = {0};
                fastq_to_fasta(&reverse_record, reverse_record_fasta, false);
                fprintf(output_forward, "%s", forward_record_fasta);
            }
            else
            {
                for (int i = 0; i < lines_to_read; i++)
                {
                    fprintf(output_forward, "%s\n", forward_record[i]);
                }
            }
            data->printed_count++;
            is_printed = true;
        }
        else
        {
            data->skipped_count++;
            is_printed = false;
        }

        // debug report (this check isn't a massive hit)
        if (debug > 1)
        {

            if (is_printed == true)
            {
                printf("Thread %d - Sequence pair %'zu PRINTED: ", thread_id, data->processed_count);
            }
            else
            {
                printf("Thread %d - Sequence pair %'zu SKIPPED: ", thread_id, data->processed_count);
            }
            printf("High (%d) count kmers: F:%d, Total unique kmers: F:%d High count ratio: F:%.2f\n",
                   depth_per_cpu,
                   seq_high_count_kmers_forward,
                   total_seq_kmers_forward,
                   high_count_ratio_fwd);
        }
        sequences_processed++;

        current_time = time(NULL);
        if (difftime(current_time, data->last_report_time) >= REPORTING_INTERVAL)
        {
            if (debug > 1)
                printf("reporting after %d seconds\n", REPORTING_INTERVAL);

            double elapsed_time = difftime(current_time, data->last_report_time);
            sequences_last_processed = data->processed_count - data->last_report_count;
            double rate = sequences_last_processed / elapsed_time;
            total_kmers = local_hash_table->used;

            float printed_improvement = (prev_printed_count == 0) ? 0 : (float)(data->printed_count - prev_printed_count) / prev_printed_count;
            float skipped_improvement = (prev_skipped_count == 0) ? 0 : (float)(data->skipped_count - prev_skipped_count) / prev_skipped_count;
            float prev_rate_improvement = (prev_rate == 0) ? 0 : (float)(rate - prev_rate) / prev_rate;
            float kmer_improvement = (prev_total_kmers == 0) ? 0 : (float)(total_kmers - prev_total_kmers) / prev_total_kmers;

            if (verbose)
                printf("Thread %d - Processing rate: %'.0f (%+.2f%%) sequences/s, processed %'zu pairs, printed: %'zu (%+.2f%%), skipped: %'zu (%+.2f%%), Unique kmers (all sequences; this thread): %'zu (%+.2f%%)\n",
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
    }

    // chunk finished processing

    if (verbose)
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

    printf("Thread %d - Processing rate: %'.0f (%+.2f%%) sequences/s, processed %'zu pairs, printed: %'zu (%+.2f%%), skipped: %'zu (%+.2f%%), Unique kmers (all sequences; this thread): %'zu (%+.2f%%)\n",
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

int multithreaded_process_files_single(thread_data_t *thread_data, mmap_file_t *forward_mmap)
{
    char stop_char = cfg.is_input_fastq == false ? '>' : '@';

    size_t max_total_kmers_in_threads = 0;
    pthread_t threads[cfg.cpus];

    // let's get start and ends for the mmapped file.
    size_t *forward_thread_starts = NULL;
    size_t *forward_thread_ends = NULL;
    forward_thread_starts = calloc(cfg.cpus, sizeof(size_t));
    forward_thread_ends = calloc(cfg.cpus, sizeof(size_t));

    if (!forward_thread_starts || !forward_thread_ends)
    {
        perror("Error allocating memory for thread positions");
        exit(EXIT_FAILURE);
    }

    // todo: bug, nb this will not work because forward and reverse can have different file sizes.
    if (cfg.cpus == 1)
    {
        if (cfg.verbose)
            printf("Single thread mode\n");

        forward_thread_ends[0] = forward_mmap->size - 1;
    }
    else
    {
        calculate_thread_positions(forward_mmap->data, forward_mmap->size, cfg.cpus, &forward_thread_starts, &forward_thread_ends);
    }

    for (int thread_number = 0; thread_number < cfg.cpus; thread_number++)
    {

        // store memory mapped file starts and stops, and data structure.
        thread_data[thread_number].forward_file_start = forward_thread_starts[thread_number];
        thread_data[thread_number].forward_file_end = forward_thread_ends[thread_number];
        thread_data[thread_number].forward_data = forward_mmap->data;

        if (!thread_data[thread_number].thread_output_forward || !thread_data[thread_number].thread_output_reverse)
        {
            fprintf(stderr, "Error opening thread-specific output files for thread %d, %s and %s\n", thread_number, thread_data[thread_number].thread_output_forward_filename, thread_data[thread_number].thread_output_reverse_filename);
            exit(EXIT_FAILURE);
        }

        // printf("forward/reverse size %zu/%zu chunk %s/%s\n", thread_data[i].forward_size, thread_data[i].reverse_size, forward_chunk_end, reverse_chunk_end);

        if (cfg.debug > 0)
            printf("Starting thread %d\n", thread_number);

        if (pthread_create(&threads[thread_number], NULL, process_thread_chunk_single, &thread_data[thread_number]) != 0)
        {

            fprintf(stderr, "Error creating thread %d\n", thread_number);
            for (int j = 0; j < thread_number; j++)
            {
                pthread_cancel(threads[j]);
                pthread_join(threads[j], NULL);
            }
            free(forward_thread_starts);
            free(forward_thread_ends);

            return 1;
        }

        sleep(1);
    }

    // close
    for (int i = 0; i < cfg.cpus; i++)
    {
        if (pthread_join(threads[i], NULL) != 0)
        {
            fprintf(stderr, "Error joining thread %d\n", i);
            free(forward_thread_starts);
            free(forward_thread_ends);

            return 1;
        }
    }

    // final report
    size_t total_processed = 0, total_printed = 0, total_skipped = 0;
    for (int thread_number = 0; thread_number < cfg.cpus; thread_number++)
    {
        total_processed += thread_data[thread_number].processed_count;
        total_printed += thread_data[thread_number].printed_count;
        total_skipped += thread_data[thread_number].skipped_count;
        max_total_kmers_in_threads = (max_total_kmers_in_threads < thread_data[thread_number].hash_table->used) ? thread_data[thread_number].hash_table->used : max_total_kmers_in_threads;

        if (cfg.doprint)
            print_kmer_table(thread_data[thread_number].hash_table, "", true, thread_number);
    }
    reporting.total_processed = total_processed;
    reporting.total_printed = total_printed;
    reporting.total_skipped = total_skipped;
    reporting.files_processed++;
    reporting.max_total_kmers = (reporting.max_total_kmers < max_total_kmers_in_threads) ? max_total_kmers_in_threads : reporting.max_total_kmers;

    printf("Cumulative file statistics: Processed %'zu, Printed %'zu, Skipped %'zu, Cumulative Max Unique Kmers in a thread: %'zu\n",
           total_processed, total_printed, total_skipped, max_total_kmers_in_threads);

    free(forward_thread_starts);
    free(forward_thread_ends);

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
    reporting.max_total_kmers = 0;
    reporting.files_processed = 0;

    memset(&cfg, 0, sizeof(struct config_t));
    if (!parse_arguments(argc, argv))
    {
        print_usage();
        return 1;
    }
    hash_table_t seed_hash = {0};

    init_hash_table(&seed_hash);
    int records_to_seed = 1 + (SEED_NUMBER / cfg.forward_file_count);

    for (int file_index = 0; file_index < cfg.forward_file_count; file_index++)
    {
        seed_kmer_hash(cfg.forward_files[file_index], records_to_seed, &seed_hash, cfg.ksize, cfg.depth);

        if (file_index < cfg.reverse_file_count)
            seed_kmer_hash(cfg.reverse_files[file_index], records_to_seed, &seed_hash, cfg.ksize, cfg.depth);
    }
    if (cfg.doprint)
        print_kmer_table(&seed_hash, "_seeds", true, -1);

    thread_data_t thread_data[cfg.cpus];

    // need a for loop for any data that will stay between runs
    for (int thread_number = 0; thread_number < cfg.cpus; thread_number++)
    {
        thread_data[thread_number].thread_id = thread_number;
        thread_data[thread_number].processed_count = 0;
        thread_data[thread_number].printed_count = 0;
        thread_data[thread_number].skipped_count = 0;
        thread_data[thread_number].last_report_time = time(NULL);
        thread_data[thread_number].last_report_count = 0;

        // if i wanted to split hash tables
        // for (int i=0;i<NUM_PARTITIONS;i++){ init_hash_table(&global_hash_tables[i]); }
        // etc

        // thread_data[thread_number].hash_table = malloc(sizeof(hash_table_t));
        // memset(thread_data[thread_number].hash_table, 0, sizeof(hash_table_t));
        // if (!thread_data[thread_number].hash_table)
        // {
        //     fprintf(stderr, "Thread %d: Memory allocation failed for hash table\n", thread_number);
        //     exit(EXIT_FAILURE);
        // }

        // init_hash_table(thread_data[thread_number].hash_table);
        thread_data[thread_number].hash_table = copy_hash_table(&seed_hash);

        // initialise the output files for each thread (otherwise we'd need to open as append)
        // every thread prints to its own file
        thread_data[thread_number].thread_output_forward_filename = NULL;
        thread_data[thread_number].thread_output_reverse_filename = NULL;

        thread_data[thread_number].thread_output_forward_filename = create_output_filename("output_forward", cfg.ksize, cfg.depth_per_cpu, thread_number, "fastq");
        thread_data[thread_number].thread_output_forward = fopen(thread_data[thread_number].thread_output_forward_filename, "w");
        if (!thread_data[thread_number].thread_output_forward)
        {
            fprintf(stderr, "Error opening file to write: %s\n", thread_data[thread_number].thread_output_forward_filename);
            exit(EXIT_FAILURE);
        }

        if (cfg.reverse_file_count != 0)
        {
            thread_data[thread_number].thread_output_reverse_filename = create_output_filename("output_reverse", cfg.ksize, cfg.depth_per_cpu, thread_number, "fastq");
            thread_data[thread_number].thread_output_reverse = fopen(thread_data[thread_number].thread_output_reverse_filename, "w");
            if (!thread_data[thread_number].thread_output_reverse)
            {
                fprintf(stderr, "Error opening file to write: %s\n", thread_data[thread_number].thread_output_reverse_filename);
                exit(EXIT_FAILURE);
            }
        }
    }
    destroy_hash_table(&seed_hash);

    // start processing files
    time_t start_time = time(NULL);

    for (int file_index = 0; file_index < cfg.forward_file_count; file_index++)
    {
        mmap_file_t forward_mmap = mmap_file(cfg.forward_files[file_index]);
        mmap_file_t reverse_mmap;
        if (file_index < cfg.reverse_file_count)
        {
            printf("Processing file pair %d of %d: %s and %s\n", file_index + 1, cfg.forward_file_count, cfg.forward_files[file_index], cfg.reverse_files[file_index]);
            reverse_mmap = mmap_file(cfg.reverse_files[file_index]);
            if (forward_mmap.data == NULL || reverse_mmap.data == NULL)
            {
                fprintf(stderr, "Error memory mapping input files\n");
                goto cleanup;
            }
        }
        else
        {
            printf("Processing single-ended file %d of %d: %s\n", file_index + 1, cfg.forward_file_count, cfg.forward_files[file_index]);
            if (forward_mmap.data == NULL)
            {
                fprintf(stderr, "Error memory mapping input files\n");
                goto cleanup;
            }
        }

        // we could skip asking user to tell us if inputs are fastq or fasta and just autodetect
        // output would be determined on whether there is at least one input fasta so would need to loop all files first and open them
        // ie not here.
        if (cfg.is_input_fastq == true && forward_mmap.data[0] != '@')
        {
            fprintf(stderr, "Input FASTQ file %s starts with %c which is not expected\n", cfg.forward_files[file_index], forward_mmap.data[0]);
            exit(EXIT_FAILURE);
        }
        else if (cfg.is_input_fastq == false && forward_mmap.data[0] != '>')
        {
            fprintf(stderr, "Input FASTA file %s starts with %c which is not expected\n", cfg.forward_files[file_index], forward_mmap.data[0]);
            exit(EXIT_FAILURE);
        }

        if (file_index < cfg.reverse_file_count)

        {
            if (cfg.is_input_fastq == true && reverse_mmap.data[0] != '@')
            {
                fprintf(stderr, "Input FASTQ file %s starts with %c which is not expected\n", cfg.reverse_files[file_index], reverse_mmap.data[0]);
                exit(EXIT_FAILURE);
            }
            else if (cfg.is_input_fastq == false && reverse_mmap.data[0] != '>')
            {
                fprintf(stderr, "Input FASTA file %s starts with %c which is not expected\n", cfg.reverse_files[file_index], reverse_mmap.data[0]);
                exit(EXIT_FAILURE);
            }
        }

        if (file_index < cfg.reverse_file_count)

        {
            if (cfg.debug)
                printf("launching paired process for file %d\n", file_index);
            // passing local copy as was having issues with pointers.
            if (multithreaded_process_files_paired(thread_data, &forward_mmap, &reverse_mmap) != 0)
            {
                fprintf(stderr, "Error processing files\n");

                munmap_file(&forward_mmap);
                munmap_file(&reverse_mmap);
                goto cleanup;
            }
            munmap_file(&forward_mmap);
            munmap_file(&reverse_mmap);
        }
        else
        {
            if (cfg.debug)
                printf("launching single process for file %d\n", file_index);
            if (multithreaded_process_files_single(thread_data, &forward_mmap) != 0)
            {
                fprintf(stderr, "Error processing files\n");

                munmap_file(&forward_mmap);
                goto cleanup;
            }
            munmap_file(&forward_mmap);
        }
    }

    // destroy_hash_table(&seed_hash);

    // closeup and free data used across all files.
    for (int thread_number = 0; thread_number < cfg.cpus; thread_number++)
    {
        fclose(thread_data[thread_number].thread_output_forward);
        free(thread_data[thread_number].thread_output_forward_filename);
        if (cfg.reverse_file_count != 0)
        {
            fclose(thread_data[thread_number].thread_output_reverse);
            free(thread_data[thread_number].thread_output_reverse_filename);
        }
        free(thread_data[thread_number].hash_table->kmer);
        free(thread_data[thread_number].hash_table);
    }

    printf("\n--- Final Report ---\n");
    printf("Processed Records: %'zu\n", reporting.total_processed);
    printf("Printed Records: %'zu\n", reporting.total_printed);
    printf("Skipped Records: %'zu\n", reporting.total_skipped);
    printf("Cumulative Max unique kmers in any thread: %'zu\n", reporting.max_total_kmers);

    // we can't get this if we have multiple threads unless we merge the tables, is it worth it?
    //    printf("Total unique kmers across all sequences: %'zu\n", reporting.total_kmers);

cleanup:

    for (int thread_number = 0; thread_number < cfg.forward_file_count; thread_number++)
    {
        free(cfg.forward_files[thread_number]);
        if (thread_number < cfg.reverse_file_count)
        {
            free(cfg.reverse_files[thread_number]);
        }
    }
    free(cfg.forward_files);
    if (cfg.reverse_file_count != 0)
        free(cfg.reverse_files);

    time_t end_time = time(NULL);
    double total_runtime = difftime(end_time, start_time);
    printf("Total runtime: %.2f seconds\n", total_runtime);
    if (reporting.total_processed > 0)
    {
        double processing_rate = reporting.total_processed / total_runtime;

        if (cfg.reverse_file_count != 0)
            printf("Overall processing rate: %'.0f sequence pairs per second\n", processing_rate);
        else
            printf("Overall processing rate: %'.0f sequences per second\n", processing_rate);
    }
    else
    {
        printf("No data processed\n");
    }
    return 0;
}
