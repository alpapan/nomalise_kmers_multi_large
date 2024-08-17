// author: alexie papanicolaou 20240814
// TODO: improvements can include considering only canonical kmers but i wrote this to process stranded RNAseq

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <getopt.h>
#include <locale.h>

#define INITIAL_CAPACITY 100000007 // near genome or transcriptome size - should be a prime number
#define MAX_SEQ_LENGTH 1024
#define MAX_FILES 1000

typedef struct
{
    uint64_t hash;
    int count;
} kmer_t;

typedef struct
{
    kmer_t *entries;
    size_t size;
    size_t capacity;
} hash_table_t;

struct config_t
{
    char *forward_files[MAX_FILES];
    char *reverse_files[MAX_FILES];
    int file_count;
    int ksize;
    int depth;
    float coverage;
    int verbose;
    char *filetype;
};

static const uint8_t base_map[256] = {
    ['A'] = 0x00, ['C'] = 0x01, ['G'] = 0x02, ['T'] = 0x03};

void init_hash_table(hash_table_t *ht)
{
    ht->entries = calloc(INITIAL_CAPACITY, sizeof(kmer_t));
    if (!ht->entries)
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    ht->size = 0;
    ht->capacity = INITIAL_CAPACITY;
}

void expand_hash_table(hash_table_t *ht)
{
    size_t new_capacity = ht->capacity * 2;
    kmer_t *new_entries = calloc(new_capacity, sizeof(kmer_t));
    if (!new_entries)
    {
        fprintf(stderr, "Memory allocation failed\n");
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
    printf("Memory expanded: %zu\n", new_capacity);
}

inline uint64_t kmer_hash(const char *seq, int k)
{
    static const uint8_t base_map[256] = {
        ['A'] = 0x00, ['C'] = 0x01, ['G'] = 0x02, ['T'] = 0x03};
    uint64_t hash = 0;
    for (int i = 0; i < k; i++)
    {
        hash = (hash << 2) | base_map[(uint8_t)seq[i]];
    }
    return hash;
}

char *create_output_filename(const char *input_filename, int k, int norm_depth)
{
    size_t len = strlen(input_filename) + 20;
    char *output_filename = malloc(len);
    if (!output_filename)
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    snprintf(output_filename, len, "%s.norm_%d_%d", input_filename, k, norm_depth);
    return output_filename;
}

inline int find_or_insert_kmer(hash_table_t *hash_table, uint64_t hash)
{
    if (hash_table->size >= hash_table->capacity * 0.75)
    {
        expand_hash_table(hash_table);
    }

    size_t index = hash % hash_table->capacity;
    size_t original_index = index;

    while (hash_table->entries[index].hash != 0 && hash_table->entries[index].hash != hash)
    {
        index = (index + 1) % hash_table->capacity;
        if (index == original_index)
        {
            // We've wrapped around the entire table
            expand_hash_table(hash_table);
            return find_or_insert_kmer(hash_table, hash);
        }
    }

    if (hash_table->entries[index].hash == 0)
    {
        hash_table->entries[index].hash = hash;
        hash_table->entries[index].count = 0;
        hash_table->size++;
    }

    return index;
}

void process_sequence(const char *seq, hash_table_t *hash_table, int K, int NORM_DEPTH, int *seq_high_count_kmers, int *total_seq_kmers)
{
    *seq_high_count_kmers = 0;
    *total_seq_kmers = 0;
    int seq_len = strlen(seq);
    for (int i = 0; i <= seq_len - K; i++)
    {
        uint64_t hash = kmer_hash(seq + i, K);
        if (hash == 0)
            continue;
        (*total_seq_kmers)++;
        int index = find_or_insert_kmer(hash_table, hash);
        if (++hash_table->entries[index].count >= NORM_DEPTH)
        {
            (*seq_high_count_kmers)++;
        }
    }
}

void print_usage(char *program_name)
{
    fprintf(stderr, "Usage: %s --forward [--forward ...] --reverse [--reverse ...] --ksize --depth [--coverage] [--verbose] [--filetype|-t fq|fa]\n", program_name);
}

int parse_arguments(int argc, char *argv[], struct config_t *cfg)
{
    cfg->coverage = 0.9;
    cfg->verbose = 0;
    cfg->file_count = 0;
    cfg->filetype = "fq";

    static struct option long_options[] = {
        {"forward", required_argument, 0, 'f'},
        {"reverse", required_argument, 0, 'r'},
        {"ksize", required_argument, 0, 'k'},
        {"depth", required_argument, 0, 'd'},
        {"coverage", required_argument, 0, 'c'},
        {"verbose", no_argument, 0, 'v'},
        {"filetype", required_argument, 0, 't'},
        {0, 0, 0, 0}};
    int opt;
    while ((opt = getopt_long(argc, argv, "F:R:K:D:C:V:T", long_options, NULL)) != -1)
    {
        switch (opt)
        {
        case 'f':
            if (cfg->file_count >= MAX_FILES)
            {
                fprintf(stderr, "Too many input files (max: %d)\n", MAX_FILES);
                return 0;
            }
            cfg->forward_files[cfg->file_count] = optarg;
            break;
        case 'r':
            if (cfg->file_count >= MAX_FILES)
            {
                fprintf(stderr, "Too many input files (max: %d)\n", MAX_FILES);
                return 0;
            }
            cfg->reverse_files[cfg->file_count] = optarg;
            cfg->file_count++;
            break;
        case 'k':
            cfg->ksize = atoi(optarg);
            break;
        case 'd':
            cfg->depth = atoi(optarg);
            break;
        case 'c':
            cfg->coverage = atof(optarg);
            break;
        case 'v':
            cfg->verbose = 1;
            break;
        case 't':
            cfg->filetype = optarg;
            break;
        default:
            return 0;
        }
    }
    if (cfg->file_count == 0 || cfg->ksize <= 0 || cfg->depth <= 0 || cfg->coverage < 0 || cfg->coverage > 1)
    {
        return 0;
    }
    return 1;
}

int main(int argc, char *argv[])
{
    setlocale(LC_ALL, "");
    struct config_t cfg = {0};
    if (!parse_arguments(argc, argv, &cfg))
    {
        print_usage(argv[0]);
        return 1;
    }

    int lines_to_read = strcmp(cfg.filetype, "fa") == 0 ? 2 : 4;
    char *output_forward_filename = create_output_filename("output_forward", cfg.ksize, cfg.depth);
    char *output_reverse_filename = create_output_filename("output_reverse", cfg.ksize, cfg.depth);
    FILE *output_forward = fopen(output_forward_filename, "w");
    FILE *output_reverse = fopen(output_reverse_filename, "w");
    if (!output_forward || !output_reverse)
    {
        fprintf(stderr, "Error opening output files\n");
        return 1;
    }

    hash_table_t hash_table;
    init_hash_table(&hash_table);
    int processed_count = 0, printed_count = 0, skipped_count = 0, prev_printed_count = 0, prev_skipped_count = 0;
    size_t total_kmers = 0;
    char forward_record[4][MAX_SEQ_LENGTH], reverse_record[4][MAX_SEQ_LENGTH];
    time_t start_time = time(NULL), last_report_time = start_time;
    int last_report_count = 0;

    for (int file_index = 0; file_index < cfg.file_count; file_index++)
    {
        FILE *forward_file = fopen(cfg.forward_files[file_index], "r");
        FILE *reverse_file = fopen(cfg.reverse_files[file_index], "r");
        if (!forward_file || !reverse_file)
        {
            fprintf(stderr, "Error opening input files\n");
            return 1;
        }

        while (1)
        {
            int forward_read = 0, reverse_read = 0;
            for (int i = 0; i < lines_to_read; i++)
            {
                if (fgets(forward_record[i], MAX_SEQ_LENGTH, forward_file) != NULL)
                {
                    forward_read++;
                }
                if (fgets(reverse_record[i], MAX_SEQ_LENGTH, reverse_file) != NULL)
                {
                    reverse_read++;
                }
            }

            if (forward_read != lines_to_read || reverse_read != lines_to_read)
            {
                break; // End of file or error
            }
            else
            {

                processed_count++;

                int forward_seq_len = strlen(forward_record[1]);
                int reverse_seq_len = strlen(reverse_record[1]);
                if (cfg.verbose)
                {
                    printf("Processing sequence pair %d from file pair %d\n", processed_count, file_index + 1);
                    printf("Forward sequence length: %d\n", forward_seq_len);
                    printf("Reverse sequence length: %d\n", reverse_seq_len);
                }
                if (forward_seq_len < cfg.ksize || reverse_seq_len < cfg.ksize)
                {
                    skipped_count++;
                    if (cfg.verbose)
                    {
                        printf("Sequence pair %d skipped due to short length\n", processed_count);
                    }
                    continue; // Skip to the next record
                }

                int seq_high_count_kmers_forward = 0, total_seq_kmers_forward = 0;
                int seq_high_count_kmers_reverse = 0, total_seq_kmers_reverse = 0;
                process_sequence(forward_record[1], &hash_table, cfg.ksize, cfg.depth, &seq_high_count_kmers_forward, &total_seq_kmers_forward);
                process_sequence(reverse_record[1], &hash_table, cfg.ksize, cfg.depth, &seq_high_count_kmers_reverse, &total_seq_kmers_reverse);
                int seq_high_count_kmers = seq_high_count_kmers_forward + seq_high_count_kmers_reverse;
                int total_seq_kmers = total_seq_kmers_forward + total_seq_kmers_reverse;
                float high_count_ratio = total_seq_kmers > 0 ? (float)seq_high_count_kmers / total_seq_kmers : 0;

                if (high_count_ratio <= cfg.coverage)
                {
                    for (int i = 0; i < lines_to_read; i++)
                    {
                        fprintf(output_forward, "%s", forward_record[i]);
                        fprintf(output_reverse, "%s", reverse_record[i]);
                    }
                    printed_count++;
                    if (cfg.verbose)
                    {
                        printf("Sequence pair %'d printed\n", processed_count);
                    }
                }
                else
                {
                    skipped_count++;
                    if (cfg.verbose)
                    {
                        printf("Sequence pair %'d skipped\n", processed_count);
                    }
                }

                if (cfg.verbose)
                {
                    // printf("Hash table size: %zu\n", hash_table.size);
                    printf("Sequence pair processed. Total kmers in pair: %d, High count kmers in pair: %d, Ratio: %.2f\n",
                           total_seq_kmers, seq_high_count_kmers, high_count_ratio);
                }

                time_t current_time = time(NULL);
                if (difftime(current_time, last_report_time) >= 60)
                {
                    double elapsed_time = difftime(current_time, last_report_time);
                    int sequences_processed = processed_count - last_report_count;
                    double rate = sequences_processed / elapsed_time;
                    printf("Processing rate: %'.0f sequences per second\n", rate);
                    last_report_time = current_time;
                    last_report_count = processed_count;
                }

                if (processed_count % 100000 == 0)
                {
                    total_kmers = hash_table.size;
                    float printed_improvement = (prev_printed_count == 0) ? 0 : (float)(printed_count - prev_printed_count) / prev_printed_count * 100;
                    float skipped_improvement = (prev_skipped_count == 0) ? 0 : (float)(skipped_count - prev_skipped_count) / prev_skipped_count * 100;

                    printf("Processed %'d sequence pairs, printed: %'d (+%.2f%%), skipped: %'d (+%.2f%%), Total kmers across all sequences: %'zu\n",
                           processed_count,
                           printed_count,
                           printed_improvement,
                           skipped_count,
                           skipped_improvement,
                           total_kmers);
                    prev_printed_count = printed_count;
                    prev_skipped_count = skipped_count;
                }
            }
        }

        fclose(forward_file);
        fclose(reverse_file);
    }

    printf("Processed Records: %d\n", processed_count);
    printf("Printed Records: %d\n", printed_count);
    printf("Skipped Records: %d\n", skipped_count);
    printf("Total kmers across all sequences: %zu\n", total_kmers);
    printf("Final hash table size: %zu\n", hash_table.size);

    time_t end_time = time(NULL);
    double total_runtime = difftime(end_time, start_time);
    printf("Total runtime: %.2f seconds\n", total_runtime);

    fclose(output_forward);
    fclose(output_reverse);
    free(hash_table.entries);
    free(output_forward_filename);
    free(output_reverse_filename);

    return 0;
}
