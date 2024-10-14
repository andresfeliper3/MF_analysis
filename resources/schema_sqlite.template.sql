-- Create organisms table
CREATE TABLE organisms (
  id INTEGER PRIMARY KEY,
  name VARCHAR,
  GCF VARCHAR,
  amount_chromosomes INT
);

-- Create chromosomes table
CREATE TABLE whole_chromosomes (
    id INTEGER PRIMARY KEY,
    name VARCHAR,
    refseq_accession_number VARCHAR,
    organism_id INTEGER REFERENCES organisms(id),
    cover_percentage REAL,
    cover REAL[],
    size INTEGER
);

CREATE TABLE region_chromosomes (
    id INTEGER PRIMARY KEY,
    name VARCHAR,
    refseq_accession_number VARCHAR,
    organism_id INTEGER REFERENCES organisms(id),
    cover_percentage REAL,
    cover REAL[],
    regions_total INTEGER,
    region_number INTEGER,
    size INTEGER,
    whole_chromosome_id INTEGER REFERENCES whole_chromosomes(id)
);

-- Create mi_grids table with foreign key constraint
CREATE TABLE whole_mi_grids (
  id INTEGER PRIMARY KEY,
  mi_grid BLOB,
  whole_chromosome_id INTEGER REFERENCES whole_chromosomes(id),
  epsilon_size REAL
);

-- Create mi_grids table with foreign key constraint
CREATE TABLE region_mi_grids (
  id INTEGER PRIMARY KEY,
  mi_grid BLOB,
  region_chromosome_id INTEGER REFERENCES region_chromosomes(id),
  epsilon_size REAL
);

-- Create chr_whole_results table with foreign key constraint
CREATE TABLE chr_whole_results (
  id INTEGER PRIMARY KEY,
  whole_chromosome_id INTEGER REFERENCES whole_chromosomes(id),
  Dq_values REAL[],
  tau_q_values REAL[],
  DDq REAL
);

-- Create chr_region_results table with foreign key constraint
CREATE TABLE chr_region_results (
  id INTEGER PRIMARY KEY,
  region_chromosome_id INTEGER REFERENCES region_chromosomes(id),
  Dq_values REAL[],
  tau_q_values REAL[],
  DDq REAL
);

-- Create repeats table for identifying general repeats
CREATE TABLE repeats (
    id INTEGER PRIMARY KEY,
    name VARCHAR,
    class_family VARCHAR,
    method_to_find_it VARCHAR
);

-- Table for identifying repeats in whole chromosomes with coordinates
CREATE TABLE recursive_repeats_whole_chromosomes (
    id INTEGER PRIMARY KEY,
    repeats_id INTEGER REFERENCES repeats(id),
    whole_chromosomes_id INTEGER REFERENCES whole_chromosomes(id),
    size INTEGER,
    largest_value INTEGER,
    coordinates VARCHAR
);


CREATE TABLE RM_repeats_whole_chromosomes (
    id INTEGER PRIMARY KEY,
    repeats_id INTEGER REFERENCES repeats(id),
    whole_chromosomes_id INTEGER REFERENCES repeats(id),
    sw_score  INTEGER,
    percentage_divergence REAL,
    percentage_deletions REAL,
    percentage_insertions REAL,
    query_begin INTEGER,
    query_end INTEGER,
    repeat_length INTEGER,
    query_left INTEGER,
    strand VARCHAR,
    repeat_begin INTEGER,
    repeat_end INTEGER,
    repeat_left INTEGER
);

CREATE TABLE gtf_genes (
    id INTEGER PRIMARY KEY,
    whole_chromosomes_id INTEGER REFERENCES whole_chromosomes(id),
    source VARCHAR,
    feature VARCHAR,
    start_position INTEGER,
    end_position INTEGER,
    length INTEGER,
    score REAL,
    strand CHAR(1),
    frame CHAR(1),
    gene_id_gtf VARCHAR,
    gene VARCHAR,
    gene_biotype VARCHAR
);

CREATE TABLE linear_repeats_whole_chromosomes (
	id INTEGER PRIMARY KEY,
    repeats_id INTEGER REFERENCES repeats(id),
    whole_chromosomes_id INTEGER REFERENCES whole_chromosomes(id),
    size INTEGER,
    count INTEGER
);

CREATE TABLE genes_containing_repeats(
	id INTEGER PRIMARY KEY,
	gtf_genes_id INTEGER REFERENCES gtf_genes(id),
	repeats_id INTEGER REFERENCES repeats(id),
	count INTEGER

);
