CREATE DATABASE IF NOT EXISTS TF_database;
USE TF_database;

CREATE TABLE IF NOT EXISTS TFlink (
    UniprotID_TF VARCHAR(20),
    UniprotID_Target VARCHAR(20),
    NCBI_GeneID_TF INT,
    NCBI_GeneID_Target INT,
    Name_TF VARCHAR(50),
    Name_Target VARCHAR(50),
    Detection_method TEXT,
    PubmedID TEXT,
    Organism VARCHAR(50),
    Source_database VARCHAR(255),
    Small_scale_evidence VARCHAR(10),
    TF_TFLink_ortho TEXT,
    TF_nonTFLink_ortho TEXT,
    Target_TFLink_ortho TEXT,
    Target_nonTFLink_ortho TEXT
);

LOAD DATA LOCAL INFILE '/home/medernoe/Ubuntu_user/Stage_M1_ubuntu/StageM1_git/Data/Transcript/Input/TFLink_Homo_sapiens_interactions_All_simpleFormat_v1.0.tsv'
INTO TABLE TFlink
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n'
IGNORE 1 LINES;
