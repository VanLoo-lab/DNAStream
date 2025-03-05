import sys 
import os
import pandas as pd 
import h5py 
import numpy as np
import json 

"""
/
 ├── SNV/                     # Shared SNV index
 │   ├── labels               #short name chr:pos:ref:alt
     |── data                 #dataframe structure containing quality scores, number of callers, etc
     |── cluster
     |── index_map             #json string for fast loading and saving
 ├── sample/                     # Shared SNV index
 │   ├── labels               #short name chr:pos:ref:alt
     |── data                 #dataframe structure containing bam file path, sample code
     |── cluster
     |── index_map
 ├── read_counts/               # Read count matrices
 │   ├── bulk/                  # Bulk sequencing read counts
 │   │   ├── variant       # SNVs x Samples (variant read counts)
 │   │   ├── total         # SNVs x Samples (total read counts)
 │   ├── lcm/                    # LCM sequencing read counts
 │   │   ├── variant       
 │   │   ├── total         
 │   ├── scdna/                  # scDNA-seq read counts
 │   │   ├── variant       
 │   │   ├── total         
 ├── metadata/                   # Metadata storage
 │   ├── sample_info              # Sample IDs
 │   ├── processing_parameters
"""



class DNAStream:
    def __init__(self, filename, verbose=True, 
            snv_dtype = np.dtype([
                    ("label", h5py.string_dtype(encoding="utf-8")),
                    ("chrom", "S5"),  # Fixed-length string (5 characters max)
                    ("pos", "i8"), 
                    ("end_pos", "i8"),    # Integer
                    ("ref_allele", "S1"),  # Single-character string
                    ("alt_allele", "S1"),  # Single-character string
                    ("hugo", "S15"), # Fixed-length string (15 characters max)
                    ("gene", "S10"), #Fixed-length string (10 characters ma)       
             ]),
            sample_dtype = np.dtype([
                    ("label", h5py.string_dtype(encoding="utf-8")),  #sample description
                    ("patient", "S10"),  # Fixed-length string (10 characters max)
                    ("source", "S10"),  # Fixed-length string (10 characters max)
                    ("location", "S15") , # Fixed-length string (15 characters max)
                    ("bam_file", h5py.string_dtype(encoding="utf-8")),  # Single-character string

            ])
    ):
        """Initialize HDF5 storage."""
        self.filename = filename
        self.verbose = verbose

        self.meta_tables = ["data", "label", "cluster", "index_map"]
        self.modalities = ["bulk", "lcm", "scdna"]
        self.indices = ["SNV", "sample"]


        self.file = h5py.File(filename, "a")  # Append mode (does not overwrite)
        

        # snv_grp = self.self.h5file.create_group("SNV")
        # snv_grp = snv_grp.attrs['label'] = np.array(shape=0, dtype=str)
        # snv_grp = snv_grp.attrs['cluster']= np.array(shape=0, dtype=np.int32)

        self.schema = {"SNV": 
                {"label" :  h5py.string_dtype(encoding="utf-8"),
                "cluster" : "i8",
                "data" : snv_dtype,
                "index_map": h5py.string_dtype("utf-8")
                },
            "sample" : 
                {
                "label":    h5py.string_dtype(encoding="utf-8") ,
                "cluster" : "i8",
                "data" : sample_dtype,
                "index_map":    h5py.string_dtype("utf-8")
                }
        }

        # Create shared SNV index if not already present
        for index in self.indices:
            for data in self.meta_tables:
                if f"{index}/{data}" not in self.file:
                    self.file.create_dataset(f"{index}/{data}", shape=(0,), maxshape=(None,), 
                                        dtype=self.schema[index][data], compression="gzip")
                    if data == "data":   
                        self.file[f"{index}/{data}"].attrs['columns'] = list(self.schema[index]["data"].names)

        # Create group structure
        for modality in self.modalities:
                group_path = f"read_counts/{modality}"
                if group_path not in self.file:
                    self.file.create_group(group_path)
                for reads in ["variant", "total"]:
                    if f"{group_path}/{reads}" not in self.file:
                        self.file.create_dataset(f"{group_path}/{reads}", shape=(0,0), maxshape=(None,None), dtype='i', compression="gzip")
                        self.file[f"{group_path}/{reads}"].dims[0].label = "SNV"
                        self.file[f"{group_path}/{reads}"].dims[1].label = "sample"


        
    

    def add_read_counts(self, fname, source):
        rc = pd.read_csv(fname)
        for _, row in rc.iterrows():
            self._add_read_count(source, row[0], row[1], row[3])

        

    
    def idx_by_label(self, label, index_name="SNV"):
        

        labs = self.file[f"{index}/label"]

        if len(labs) ==0:
            return None

        indices = np.where(labs[:] ==label)[0]
        if len(indices) ==0:
            return None 
        elif len(indices) == 1:
            return indices[0]
        else:
            raise ValueError("SNV label is associated with multiple indices and could not be added.")
     
    

    def label_by_idx(self, idx, index_name="SNV"):

        return self.file[f"{index_name}/label"][idx]


    def _resize_all(self, m=None,n=None ):

        if m:
            for snv_data in self.meta_tables:
                group_path = f"SNV/{snv_data}"
            
                self.file[group_path].resize((m,))
        else:
            m= self.file[f"SNV/label"].shape[0]

        if n:
            for sample_data in self.meta_tables:
                group_path = [f"sample/{sample_data}"]
                self.file[group_path].resize((m,))
        else:
            n = self.file[f"sample/label"].shape[0]


        for modality in self.modalities:
            group_path = f"read_counts/{modality}"
            for reads in ["variant", "total"]:
                
                mat = self.file[f"{group_path}/{reads}"]
                mat.resize((m, n))


    def add_snv(self, label, cluster=None, data=None, overwrite=False):
        return self.add_item(label, index="SNV", cluster=cluster, data=data, overwrite=overwrite)
    

    def add_sample(self, label, cluster=None, data=None, overwrite=False):
        return self.add_item(label, index="sample", cluster=cluster, data=data, overwrite=overwrite)


    def add_item(self, label:str, index_name, index_dict=None, cluster=None, data=None, overwrite=False):
        label = label.encode("utf-8")
        idx = self.idx_by_label(label, index_name, index_dict)

        if not idx:
            
            new_idx = self.file[f"{index}/label"].shape[0]
            if index == "SNV":
                self._resize_all(m =new_idx+1)
            else:
                self._resize_all(n=new_idx+1)
        else:
            if overwrite:
                new_idx = idx 
            else:
                print(f"{label} exists in {index} index, use overwrite=True to overwrite metadata.")
                return idx 
        
        self.file[f"{index}/label"][new_idx] = label
        self.file[f"{index}/index_map"][label] = new_idx
        if data:
            self.file[f"{index}/data"][new_idx] = data 
        
        if cluster:
            self.file[f"{index}/cluster"][new_idx] = cluster 


        return new_idx


    def load_snv_index(self):
        return self._load_index("SNV")
    
    def load_sample_index(self):
        return self._load_index("sample")
    

    def save_snv_index(self, index_dict):
        return self._save_index(index_dict, "SNV")
    
    def save_sample_index(self, index_dict):
        return self._save_index(index_dict, "sample")
    
      

    def _load_index(self, index_name):
        """Load the label -> index mapping into memory."""
        if f"{index_name}/index_map" in self.file:
            index_data = self.file[f"{index_name}/index_map"][()]
            return json.loads(index_data[0]) if len(index_data) > 0 else {}
        return {}
    

    def _save_index(self, index_dict, index_name):
        """Save the label → index mapping into HDF5."""
        index_json = json.dumps(index_dict)  # Convert dictionary to JSON string
        self.file[f"{index_name}/index_map"].resize((1,))  # Ensure space in dataset
        self.file[f"{index_name}/index_map"][0] = index_json  # Store JSON string in HDF5

    def __add_read_count(self, source, snv_label, sample_label, var=0, total=0):
        snv_idx = self.add_snv(snv_label)
        sample_idx = self.add_sample(sample_label)
    
        dat = self.file[f"read_counts/{source}"]
        dat['variant'][snv_idx,sample_idx] = var 
        dat['total'][snv_idx,sample_idx] = total


    def batch_add_snvs(self, labels, index_dict):
      
        indices = []
 
        for lab in labels:
            if lab in index_dict:
                indices.append(index_dict[lab])
            else:
                next_idx = len(index_dict)
                indices.append(next_idx)
                index_dict[lab] = next_idx

        
        if self.file[f"SNV/label"].shape[0 ]!= len(index_dict):
            self._resize_all(m=len(index_dict))
    

        self.file[f"SNV/label"][indices] = labels

        return indices
    
    def get_snv_data(self, indices=None):

        return self._get_data(dataset_name="SNV/data", indices=indices)
    


        
    def get_sample_data(self, indices=None):

        return self._get_data(dataset_name="sample/data", indices=indices)



    
    def _get_data(self,  dataset_name, indices=None):
        dataset = self.file[dataset_name]
        columns = dataset.attrs['columns']

        if indices:
            dataset = dataset[indices, :]
        else:
            dataset = dataset[:]
        
   
        
  
        df = pd.DataFrame(dataset, columns=columns)
        for col in df.select_dtypes(include=["object"]):  # Select only object (bytes) columns
            if df[col].dtype == object:  # Pandas converts 'SXX' to 'object'
                df[col] = df[col].apply(lambda x: x.decode("utf-8") if isinstance(x, bytes) else x)
    
        return df




    def add_snv_data(self, indices, df):
        self._add_data(indices, df, dataset_name="SNV/data")
    

    def add_sample_data(self, indices, df):
        self._add_data(indices, df, dataset_name="sample/data")

    
    def _add_data(self, indices, df, dataset_name):
            # Resize dataset if needed
            dataset = self.file[dataset_name]
            max_index = max(indices) + 1
            if dataset.shape[0] < max_index:
                self.file.close()
                raise IndexError("Invalid indices passed to add_data method, check indices and try again")

            # Convert DataFrame to structured NumPy array
            structured_data = np.zeros(len(indices), dtype=dataset.dtype)
            for col in df.columns:
                structured_data[col] = df[col].to_numpy()

            # Update data at provided indices
            dataset[indices] = structured_data
                
    def update_snvs_from_maf(self, fname, 
                                 missing_values= ["Unknown", "Na", "N/A", "na", "nan", 
                                                    "NaN", "NAN", "NONE", "None", "", "__UNKNOWN__"],
                                 required_cols =["Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
                                             "Reference_Allele", "Tumor_Seq_Allele2", "Entrez_Gene_Id"] ):
        """
        Create an SNV index by reading in the MAF file

        fname: str path to MAF file
        """  

        #read MAF file and extract key info
        maf = pd.read_table(fname, low_memory=False) 
      
        missing_cols = set(required_cols) - set(maf.columns)

        maf.replace(missing_values, pd.NA, inplace=True)
        if missing_cols:
            maf = maf.reindex(columns=maf.columns.tolist() + list(missing_cols), fill_value=pd.NA)
        
        
        
        maf["label"] = maf.iloc[:, :4].astype(str).agg(":".join, axis=1)

        column_dict = {
            "label" : "label",
            "Chromosome" : "chrom",
            "Start_Position" : "pos",
            "End_Position" : "end_pos",
            "Reference_Allele" : "ref_allele",
            "Tumor_Seq_Allele2" : "alt_allele",
            "Hugo_Symbol" : "hugo",
            "Entrez_Gene_Id": "gene"
        }

        maf.rename(columns = column_dict, inplace=True)

        snv_data = maf[[val for _, val in column_dict.items()]]
        if self.verbose:
            print(snv_data.head())

    
        snv_labels = maf["label"].tolist() 

        

        try: 

                snv_index = self.load_snv_index()
        
                #returns the snv indices of snv_labels
                indices = self.batch_add_snvs(snv_labels, snv_index)
        
                self.add_snv_data(indices, snv_data)
                
                #update the index in the HDF5 file
                self.save_snv_index(snv_index)

        except:
       
            self.close()
            raise ValueError("Error adding batch SNVs")

        if self.verbose:
            print(f"{len(indices)} new SNVs added to index")


    

    
    
    def close(self):
        self.file.close()


ds = DNAStream(filename="temp.h5", verbose=True)
sample = "GEM2.14_LR_1"
maf_file = f"data_summary/WGS_MUTATION/Somatic/2outof3_SNV/{sample}_SNVs_2outof3.maf"
indices = ds.update_snvs_from_maf(maf_file)
print(ds.get_snv_data().head())
ds.close()




             
                  
             

    # def loadReadCounts(self,  filedict: dict):
    #     if len(self.ReadCounts) ==0:
    #         self.SNVIdx = Index()
    #         self.SampleIdx = Index()
    #     for key,val in filedict.items():
    #         rc = ReadCounts()
    #         rc.load(val)
    #         if key in self.ReadCounts:
    #             print(f"Warning: overwriting read counts layer {key}!")

    #         self.ReadCounts[key] = rc

    # def RankTrees(self):
    #     pass 






# class Index:
#     def __init__(self, idx=None):
        
    
#         if idx is None:
#             self.idx = {}
#             self.reverse_idx = {}
#         else:
#             self.reverse_idx = {val: key for key,val in idx.items()}
        
#         self.max_idx = max(self.reverse_idx.keys())
        
#     def add(self, val):
#         if val not in self.idx:
#             self.max_idx += 1



#     def new_idx(self):
#         max()


# class ReadCounts:
#     def __init__(self):

#         self.variant = None 
#         self.total = None
#         self.SNVIdx = None 

        

    # def load(fname, snv_to_idx:Index = None, sample_to_idx: Index = None  ):
        

    #     dat  = pd.read_csv(fname)

        






        