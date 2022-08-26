#!/usr/bin/env python

import pandas as pd
import math

# Load technical test data from excel file stored on hard-disk. Expect the script to be same folder as data.
xls = pd.ExcelFile('Technical Test - Data Wrangling.xlsx')
df_patient_clinical_data = pd.read_excel(xls, 'Patient_clinical_data')
df_tissue_sample_metadata = pd.read_excel(xls, 'Tissue Sample Metadata')
df_serum_protein_data = pd.read_excel(xls, 'Serum Protein data')
df_rna_seq_rpkm = pd.read_excel(xls, 'RNA-seq (RPKM)')


# Data Wrangling


# Data cleaning for patient clinical data
# Step 1: M, F to be converted to Male, Female in patient data
df_patient_clinical_data['Sex']. replace('M', 'MALE', inplace=True)
df_patient_clinical_data['Sex']. replace('F', 'FEMALE', inplace=True)

# Step 2: Round off age
df_patient_clinical_data['Age'] = df_patient_clinical_data['Age'].round(0).astype(int)

# Step 3: Create column Unique Patient ID and renaming Patient Number column
df_patient_clinical_data["Unique_Patient_ID"] = df_patient_clinical_data["Study_ID"].astype(str) +"_"+ df_patient_clinical_data["Patient  Number"].astype(str)
df_patient_clinical_data.rename(columns = {'Patient  Number':'Patient_ID'}, inplace = True)



# Data cleaning for tissue sample metadata
# Step 1: Map the sample type
df_tissue_sample_metadata["Sample type"].replace({"Normal": "NORMAL", "Liver Tumor": "PRIMARY", "Metastic Lung": "METASTATIC"}, inplace=True)

# Step 3: Create column Unique Patient ID and renaming Patient Number column
df_tissue_sample_metadata.rename(columns = {'Patient  Number':'Patient_ID'}, inplace = True)



# Data cleaning for Serum Protein data
# Step 1: Convert the results for IL6 and IL6R to float and making the non-numeric to NaN
df_serum_protein_data["Serum IL-6 Receptor (mg/L)"] = pd.to_numeric(df_serum_protein_data["Serum IL-6 Receptor (mg/L)"], errors = 'coerce')
df_serum_protein_data["Serum IL-6 (g/L)"] = pd.to_numeric(df_serum_protein_data["Serum IL-6 (g/L)"], errors = 'coerce')

# Step 2: Make the unit of IL6R from mg/L to g/L
df_serum_protein_data["Serum IL-6 Receptor (mg/L)"] = df_serum_protein_data["Serum IL-6 Receptor (mg/L)"]/1000
df_serum_protein_data["Units"] = "g/L"

# Step 3: Rename the columns for IL6 and IL6R measurements
df_serum_protein_data.rename(columns = {'Serum IL-6 (g/L)':'IL6', 'Serum IL-6 Receptor (mg/L)':'IL6R'}, inplace = True)




# Data generation
# Generate data for tissue sample data

# Merging df_patient_clinical_data with df_tissue_sample_metadata 
result_patient_tissue = pd.merge(df_patient_clinical_data, df_tissue_sample_metadata, on=['Patient_ID'])

# Creating dataframe for final report
report_df = pd.DataFrame(columns = ['Study_ID', 'Patient_ID', 'Unique_Patient_ID','Sex', 'Age',
       'Sample_ID', 'Sample_General_Patology', 'Material_type', 'Gene_Symbol', 'Result', 'Result_Units', 'Status'])

# Iterating the tissue sample data dataframe and find the result type of gene type based on sample_id
for index, row in result_patient_tissue.iterrows():
    result_patient_tissue_entry = {'Study_ID': row ['Study_ID'], 'Patient_ID' : row ['Patient_ID'],
    'Unique_Patient_ID': row['Unique_Patient_ID'],'Sex': row['Sex'], 'Age' : row['Age'],'Sample_ID' : row['Sample'],
    'Sample_General_Patology' : row ['Sample type'], 'Material_type':  row['Material'], 'Result_Units' : 'RPKM'}
    
    sample_id = row['Sample']
    result_patient_tissue_gene_entry = {}
    result_patient_tissue_gene_entry.update(result_patient_tissue_entry)
    if sample_id in df_rna_seq_rpkm.columns:
        gene_values = df_rna_seq_rpkm.loc[:, ['GeneID', sample_id]]
        for index, value in gene_values.iterrows():
            result_patient_tissue_gene_entry['Gene_Symbol'] = value['GeneID']
            result_patient_tissue_gene_entry['Result'] = value[sample_id]
            # The non-number is considered to be not-done status, 0 is considered to be expected status, hence NA.
            result_patient_tissue_gene_entry['Status'] = 'NOT DONE' if math.isnan(value[sample_id]) else 'NA'
            report_df = report_df.append(result_patient_tissue_gene_entry, ignore_index=True)
    else:
        # Since the data is not existant, hence Gene Id is NA	
        result_patient_tissue_gene_entry['Gene_Symbol'] = 'NA'
        result_patient_tissue_gene_entry['Result'] = 'NULL'
        # Since the data is not present here, hence status is NOT DONE.
        result_patient_tissue_gene_entry['Status'] = 'NOT DONE'
        report_df = report_df.append(result_patient_tissue_gene_entry, ignore_index=True)




# Generate data for serum sample data
# Convert the result data stored for IL6 and IL6R in two columns to separate rows
report_serum_df = pd.DataFrame(columns = ['Patient_ID', 'Sample_ID', 'Gene_Symbol', 'Result', 'Status'])
# iterate on serum dataframe
for index, row in df_serum_protein_data.iterrows():
    Il6_dict = {'Patient_ID': row['Patient'],
                'Sample_ID' : row['Sample'],
                'Gene_Symbol' : 'IL6',
                'Result' : 'NULL' if math.isnan(row['IL6']) else row['IL6'],
                'Status' : 'NOT DONE' if math.isnan(row['IL6']) else 'NA',
                'Result_Units' : row['Units']
               }
    Il6R_dict = {'Patient_ID': row['Patient'],
                'Sample_ID' : row['Sample'],
                'Gene_Symbol' : 'IL6R',
                'Result' : 'NULL' if math.isnan(row['IL6R']) else row['IL6R'],
                'Status' : 'NOT DONE' if math.isnan(row['IL6R']) else 'NA',
                'Result_Units' : row['Units']
               }
    report_serum_df = report_serum_df.append(Il6_dict, ignore_index=True)
    report_serum_df = report_serum_df.append(Il6R_dict, ignore_index=True)

# Adding constant column names 
report_serum_df['Material_type'] = 'SERUM'
report_serum_df['Sample_General_Patology'] = 'NA'

# Merge patient data to serum data
report_serum_patient_clinial_df = pd.merge(df_patient_clinical_data, report_serum_df, on=['Patient_ID'])



# Append the serum data to already generated tissue data in report dataframe
report_df = report_df.append(report_serum_patient_clinial_df)



# Save the report data frame as CSV in same folder as script and data
report_df.to_csv('report_technical_test.csv', index=False)




