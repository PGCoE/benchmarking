#! /usr/env/Rscript

library("ape")
library("readr")

# obtain arguments
args <- commandArgs(trailingOnly=TRUE)

# read the input file
in_data <- read.csv(args[1])

# create an ID column derived from the names of the pairs
in_data$ID <- paste(in_data$seq1, "|", in_data$seq2)

# extract columns for PCoA
extracted_cols <- c("match", "match_iupac", "ins", "del", "sub", "invalid", 
                    "valid", "del_terminal", "del_internal", "ins_terminal", 
                    "ins_internal", "ID")
pcoa_input <- in_data[, extracted_cols]
rownames(pcoa_input) <- pcoa_input$ID
pcoa_input$ID <- NULL

# create the distance matrix
dist_matrix <- dist(pcoa_input, method="euclidean")

# run PCoA and extract coordinates
pcoa_output <- pcoa(dist_matrix, correction="none", rn=rownames(dist_matrix))
pcoa_vec1 <- pcoa_output$vectors[, 1]
pcoa_vec2 <- pcoa_output$vectors[, 2]

# bind coordinates to input dataframe
pcoa_input$pco1 <- pcoa_vec1
pcoa_input$pco2 <- pcoa_vec2

# run correlation analysis between input variables and PCoA axes
# (my lack of R skills is showing)
cor_df <- data.frame(matrix(ncol=3, nrow=0))
colnames(cor_df) <- c("var", "pco1_r2", "pco2_r2")

lm_match_1 <- lm(match ~ pco1, data=pcoa_input)
lm_match_2 <- lm(match ~ pco2, data=pcoa_input)
r2_match_1 <- summary(lm_match_1)$r.squared
r2_match_2 <- summary(lm_match_2)$r.squared
new_row <- data.frame(
    var="match", pco1_r2=r2_match_1, pco2_r2=r2_match_2 
)
cor_df <- rbind(cor_df, new_row)

lm_ins_1 <- lm(ins ~ pco1, data=pcoa_input)
lm_ins_2 <- lm(ins ~ pco2, data=pcoa_input)
r2_ins_1 <- summary(lm_ins_1)$r.squared
r2_ins_2 <- summary(lm_ins_2)$r.squared
new_row <- data.frame(
    var="ins", pco1_r2=r2_ins_1, pco2_r2=r2_ins_2 
)
cor_df <- rbind(cor_df, new_row)

lm_del_1 <- lm(del ~ pco1, data=pcoa_input)
lm_del_2 <- lm(del ~ pco2, data=pcoa_input)
r2_del_1 <- summary(lm_del_1)$r.squared
r2_del_2 <- summary(lm_del_2)$r.squared
new_row <- data.frame(
    var="del", pco1_r2=r2_del_1, pco2_r2=r2_del_2 
)
cor_df <- rbind(cor_df, new_row)

lm_sub_1 <- lm(sub ~ pco1, data=pcoa_input)
lm_sub_2 <- lm(sub ~ pco2, data=pcoa_input)
r2_sub_1 <- summary(lm_sub_1)$r.squared
r2_sub_2 <- summary(lm_sub_2)$r.squared
new_row <- data.frame(
    var="sub", pco1_r2=r2_sub_1, pco2_r2=r2_sub_2 
)
cor_df <- rbind(cor_df, new_row)

lm_invalid_1 <- lm(invalid ~ pco1, data=pcoa_input)
lm_invalid_2 <- lm(invalid ~ pco2, data=pcoa_input)
r2_invalid_1 <- summary(lm_invalid_1)$r.squared
r2_invalid_2 <- summary(lm_invalid_2)$r.squared
new_row <- data.frame(
    var="invalid", pco1_r2=r2_invalid_1, pco2_r2=r2_invalid_2 
)
cor_df <- rbind(cor_df, new_row)

lm_match_iupac_1 <- lm(match_iupac ~ pco1, data=pcoa_input)
lm_match_iupac_2 <- lm(match_iupac ~ pco2, data=pcoa_input)
r2_match_iupac_1 <- summary(lm_match_iupac_1)$r.squared
r2_match_iupac_2 <- summary(lm_match_iupac_2)$r.squared
new_row <- data.frame(
    var="match_iupac", pco1_r2=r2_match_iupac_1, pco2_r2=r2_match_iupac_2 
)
cor_df <- rbind(cor_df, new_row)

lm_valid_1 <- lm(valid ~ pco1, data=pcoa_input)
lm_valid_2 <- lm(valid ~ pco2, data=pcoa_input)
r2_valid_1 <- summary(lm_valid_1)$r.squared
r2_valid_2 <- summary(lm_valid_2)$r.squared
new_row <- data.frame(
    var="valid", pco1_r2=r2_valid_1, pco2_r2=r2_valid_2 
)
cor_df <- rbind(cor_df, new_row)

lm_del_terminal_1 <- lm(del_terminal ~ pco1, data=pcoa_input)
lm_del_terminal_2 <- lm(del_terminal ~ pco2, data=pcoa_input)
r2_del_terminal_1 <- summary(lm_del_terminal_1)$r.squared
r2_del_terminal_2 <- summary(lm_del_terminal_2)$r.squared
new_row <- data.frame(
    var="del_terminal", pco1_r2=r2_del_terminal_1, pco2_r2=r2_del_terminal_2 
)
cor_df <- rbind(cor_df, new_row)

lm_ins_terminal_1 <- lm(ins_terminal ~ pco1, data=pcoa_input)
lm_ins_terminal_2 <- lm(ins_terminal ~ pco2, data=pcoa_input)
r2_ins_terminal_1 <- summary(lm_ins_terminal_1)$r.squared
r2_ins_terminal_2 <- summary(lm_ins_terminal_2)$r.squared
new_row <- data.frame(
    var="ins_terminal", pco1_r2=r2_ins_terminal_1, pco2_r2=r2_ins_terminal_2 
)
cor_df <- rbind(cor_df, new_row)

lm_del_internal_1 <- lm(del_internal ~ pco1, data=pcoa_input)
lm_del_internal_2 <- lm(del_internal ~ pco2, data=pcoa_input)
r2_del_internal_1 <- summary(lm_del_internal_1)$r.squared
r2_del_internal_2 <- summary(lm_del_internal_2)$r.squared
new_row <- data.frame(
    var="del_internal", pco1_r2=r2_del_internal_1, pco2_r2=r2_del_internal_2 
)
cor_df <- rbind(cor_df, new_row)

lm_ins_internal_1 <- lm(ins_internal ~ pco1, data=pcoa_input)
lm_ins_internal_2 <- lm(ins_internal ~ pco2, data=pcoa_input)
r2_ins_internal_1 <- summary(lm_ins_internal_1)$r.squared
r2_ins_internal_2 <- summary(lm_ins_internal_2)$r.squared
new_row <- data.frame(
    var="ins_internal", pco1_r2=r2_ins_internal_1, pco2_r2=r2_ins_internal_2 
)
cor_df <- rbind(cor_df, new_row)

cat(format_csv(cor_df))