-- Create organisms table
CREATE TABLE organisms (
  id INT AUTO_INCREMENT PRIMARY KEY,
  name VARCHAR(255),
  GCF VARCHAR(255),
  amount_chromosomes INT
);

-- Create whole_chromosomes table
CREATE TABLE whole_chromosomes (
    id INT AUTO_INCREMENT PRIMARY KEY,
    name VARCHAR(255),
    refseq_accession_number VARCHAR(255),
    organism_id INT,
    cover_percentage FLOAT,
    cover FLOAT,
    size INT,
    FOREIGN KEY (organism_id) REFERENCES organisms(id)
);

-- Create region_chromosomes table
CREATE TABLE region_chromosomes (
    id INT AUTO_INCREMENT PRIMARY KEY,
    name VARCHAR(255),
    refseq_accession_number VARCHAR(255),
    organism_id INT,
    cover_percentage FLOAT,
    cover FLOAT,
    regions_total INT,
    region_number INT,
    size INT,
    whole_chromosome_id INT,
    FOREIGN KEY (organism_id) REFERENCES organisms(id),
    FOREIGN KEY (whole_chromosome_id) REFERENCES whole_chromosomes(id)
);

-- Create mi_grids table with foreign key constraint
CREATE TABLE whole_mi_grids (
  id INT AUTO_INCREMENT PRIMARY KEY,
  mi_grid LONGBLOB,
  whole_chromosome_id INT,
  epsilon_size FLOAT,
  FOREIGN KEY (whole_chromosome_id) REFERENCES whole_chromosomes(id)
);

-- Create mi_grids table with foreign key constraint
CREATE TABLE region_mi_grids (
  id INT AUTO_INCREMENT PRIMARY KEY,
  mi_grid LONGBLOB,
  region_chromosome_id INT,
  epsilon_size FLOAT,
  FOREIGN KEY (region_chromosome_id) REFERENCES region_chromosomes(id)
);


-- Create chr_whole_results table with foreign key constraint
CREATE TABLE chr_whole_results (
  id INT AUTO_INCREMENT PRIMARY KEY,
  whole_chromosome_id INT,
  Dq_values FLOAT,
  tau_q_values FLOAT,
  DDq FLOAT,
  FOREIGN KEY (whole_chromosome_id) REFERENCES whole_chromosomes(id)
);

-- Create chr_region_results table with foreign key constraint
CREATE TABLE chr_region_results (
  id INT AUTO_INCREMENT PRIMARY KEY,
  region_chromosome_id INT,
  Dq_values FLOAT,
  tau_q_values FLOAT,
  DDq FLOAT,
  FOREIGN KEY (region_chromosome_id) REFERENCES region_chromosomes(id)
);

-- Create repeats table for identifying general repeats
CREATE TABLE repeats (
    id INT AUTO_INCREMENT PRIMARY KEY,
    name VARCHAR(255),
    class_family VARCHAR(255),
    method_to_find_it VARCHAR(255)
);

-- Table for identifying repeats in whole chromosomes with coordinates
CREATE TABLE recursive_repeats_whole_chromosomes (
    id INT AUTO_INCREMENT PRIMARY KEY,
    repeats_id INT,
    whole_chromosomes_id INT,
    size INT,
    largest_value INT,
    coordinates VARCHAR(255),
    FOREIGN KEY (repeats_id) REFERENCES repeats(id),
    FOREIGN KEY (whole_chromosomes_id) REFERENCES whole_chromosomes(id)
);

-- Create RM_repeats_whole_chromosomes table with foreign key constraint
CREATE TABLE RM_repeats_whole_chromosomes (
    id INT AUTO_INCREMENT PRIMARY KEY,
    repeats_id INT,
    whole_chromosomes_id INT,
    sw_score INT,
    percentage_divergence FLOAT,
    percentage_deletions FLOAT,
    percentage_insertions FLOAT,
    query_begin INT,
    query_end INT,
    repeat_length INT,
    query_left INT,
    strand VARCHAR(255),
    repeat_begin INT,
    repeat_end INT,
    repeat_left INT,
    FOREIGN KEY (repeats_id) REFERENCES repeats(id),
    FOREIGN KEY (whole_chromosomes_id) REFERENCES whole_chromosomes(id)
);

-- Create gtf_genes table
CREATE TABLE gtf_genes (
    id INT AUTO_INCREMENT PRIMARY KEY,
    whole_chromosomes_id INT,
    source VARCHAR(255),
    feature VARCHAR(255),
    start_position INT,
    end_position INT,
    length INT,
    score FLOAT,
    strand CHAR(1),
    frame CHAR(1),
    gene_id_gtf VARCHAR(255),
    gene VARCHAR(255),
    gene_biotype VARCHAR(255),
    category
    FOREIGN KEY (whole_chromosomes_id) REFERENCES whole_chromosomes(id)
);


CREATE TABLE linear_repeats_whole_chromosomes (
    id INT PRIMARY KEY,
    repeats_id INT,
    whole_chromosomes_id INT,
    size BIGINT,
    count INT,
    FOREIGN KEY (repeats_id) REFERENCES repeats(id),
    FOREIGN KEY (whole_chromosomes_id) REFERENCES whole_chromosomes(id)
);

CREATE TABLE linear_repeats_region_chromosomes (
    id INT PRIMARY KEY AUTO_INCREMENT,
    repeats_id INT,
    region_chromosomes_id INT,
    size INT,
    count INT,
    FOREIGN KEY (repeats_id) REFERENCES repeats(id),
    FOREIGN KEY (region_chromosomes_id) REFERENCES region_chromosomes(id)
);

CREATE TABLE genes_containing_repeats (
    id INT AUTO_INCREMENT PRIMARY KEY,
    gtf_genes_id INT,
    repeats_id INT,
    count INT,
    FOREIGN KEY (gtf_genes_id) REFERENCES gtf_genes(id),
    FOREIGN KEY (repeats_id) REFERENCES repeats(id)
);