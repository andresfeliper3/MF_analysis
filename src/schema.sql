-- Create organisms table
CREATE TABLE organisms (
  id INTEGER PRIMARY KEY,
  name VARCHAR,
  GCF VARCHAR,
  amount_chromosomes INT
);

-- Create chromosomes table
CREATE TABLE chromosomes (
    id INTEGER PRIMARY KEY,
    name VARCHAR,
    organism_id INTEGER REFERENCES organisms(id),
    cover_percentage REAL,
    cover REAL[]
)

-- Create mi_grids table with foreign key constraint
CREATE TABLE mi_grids (
  id INTEGER PRIMARY KEY,
  mi_grid BLOB,
  chromosome_id INTEGER REFERENCES chromosomes(id),
  epsilon_size REAL
);

-- Create chr_whole_results table with foreign key constraint
CREATE TABLE chr_whole_results (
  id INTEGER PRIMARY KEY,
  chromosome_id INTEGER REFERENCES chromosomes(id),
  Dq_values REAL[],
  tau_q_values REAL[],
  DDq REAL
);

-- Create chr_region_results table with foreign key constraint
CREATE TABLE chr_region_results (
  id INTEGER PRIMARY KEY,
  regions_number INTEGER,
  chromosome_id INTEGER REFERENCES chromosomes(id),
  Dq_values REAL[],
  tau_q_values REAL[],
  DDq REAL
);
