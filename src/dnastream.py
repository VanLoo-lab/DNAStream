import sys 
import os
import getpass
import socket
import pathlib
import time
from datetime import datetime
import functools
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
     |-- log
 ├── sample/                     # Shared SNV index
 │   ├── labels               #short name chr:pos:ref:alt
     |── data                 #dataframe structure containing bam file path, sample code
     |── cluster
     |── index_map
     |-- log
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




def timeit(func):
    """Decorator to measure execution time of a function."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.perf_counter()  # Start timer
        result = func(*args, **kwargs)    # Run the function
        end_time = time.perf_counter()    # End timer
        elapsed_time = end_time - start_time
        print(f"⏱ Function '{func.__name__}' took {elapsed_time:.4f} seconds")
        return result
    return wrapper


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

        log_dtype = np.dtype([
            ("timestamp", h5py.string_dtype(encoding="utf-8")),
            ("number", "i8"),
            ("index_size_before", "i8"),
            ("index_size_after", "i8"),
            ("operation", "S15"),
            ("user", "S15"),
            ("hostname", "S15"),
            ("file",  h5py.string_dtype(encoding="utf-8"))
        ])

        self.meta_tables = ["data", "label", "cluster", "index_map", "log"]
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
                "index_map": h5py.string_dtype("utf-8"),
                "log" : log_dtype,
                },
            "sample" : 
                {
                "label":    h5py.string_dtype(encoding="utf-8") ,
                "cluster" : "i8",
                "data" : sample_dtype,
                "index_map":    h5py.string_dtype("utf-8"),
                "log" : log_dtype
                }
        }

        # Create shared SNV index if not already present
        for index in self.indices:
            for data in self.meta_tables:
                if f"{index}/{data}" not in self.file:
                    self.file.create_dataset(f"{index}/{data}", shape=(0,), maxshape=(None,), 
                                        dtype=self.schema[index][data], compression="gzip")
                    if data in ["data", "log"]:   
                        self.file[f"{index}/{data}"].attrs['columns'] = list(self.schema[index][data].names)

        # Create group structure
        for modality in self.modalities:
                group_path = f"read_counts/{modality}"
                if group_path not in self.file:
                    self.file.create_group(group_path)
                for reads in ["variant", "total"]:
                    if f"{group_path}/{reads}" not in self.file:
                        self.file.create_dataset(f"{group_path}/{reads}", shape=(0,0), maxshape=(None,None), dtype='i', 
                                                 compression="gzip", chunks=(1, 5000))
                        self.file[f"{group_path}/{reads}"].dims[0].label = "SNV"
                        self.file[f"{group_path}/{reads}"].dims[1].label = "sample"


        
    def __str__(self):
        m = self.file["SNV/label"].shape[0]
        n = self.file["sample/label"].shape[0]

        mystr = f"DNAStream object with {m} SNVs and {n} samples" 
        mystr += f"\nHDF5 File: {self.filename}"
      
        return mystr
    
    @timeit
    def add_read_counts(self, fname, source, location=None):
        try:
            rc = pd.read_csv(fname, names=["snv","sample", "var", "total"], header=None, skiprows=1)
       

            samples = rc["sample"].unique().tolist()

            snvs = rc["snv"].unique().tolist()

            snv_idx = self.batch_add_snvs(snvs, source_file=fname)

            sample_idx = self.batch_add_samples(samples, source_file=fname)
            sample_indices = list(sample_idx.values())
            if source:
                self._update_value(sample_indices, "sample/data", "source", source)
            
            if location:
                self._update_value(sample_indices, "sample/data", "location", location)

            

                
            # Map indices for SNVs and samples
            rc["snv_idx"] = rc["snv"].map(snv_idx)
            rc["sample_idx"] = rc["sample"].map(sample_idx)

            # Convert to NumPy arrays for efficient updates
            snv_indices_arr = rc["snv_idx"].to_numpy()
            sample_indices_arr = rc["sample_idx"].to_numpy()
            var_counts = rc["var"].to_numpy()
            total_counts = rc["total"].to_numpy()

            # var = self.file[f"read_counts/{source}/variant"][:]
            # total = self.file[f"read_counts/{source}/total"][:]


            # var[np.ix_(snv_indices_arr, sample_indices_arr)] = var_counts
            # total[np.ix_(snv_indices_arr, sample_indices_arr)] = total_counts

            # self.file[f"read_counts/{source}/variant"] = var 
            # self.file[f"read_counts/{source}/total"] = total

            

            # print(self)
               # **Sort indices to satisfy HDF5 fancy indexing rules**
            sorted_order = np.lexsort((sample_indices_arr, snv_indices_arr))  
            snv_indices_arr = snv_indices_arr[sorted_order]
            sample_indices_arr = sample_indices_arr[sorted_order]
            var_counts = var_counts[sorted_order]
            total_counts = total_counts[sorted_order]


            # Update dataset in batch instead of looping
            dat = self.file[f"read_counts/{source}"]
            unique_snv_indices = np.unique(snv_indices_arr)  # Get unique SNV indices
            
    

            for snv in unique_snv_indices:
                mask = snv_indices_arr == snv  # Select all entries for this snv index
                dat["variant"][snv, sample_indices_arr[mask]] = var_counts[mask]
                dat["total"][snv, sample_indices_arr[mask]] = total_counts[mask]




        except Exception as e:
            self.close()
            raise Exception(e)


        
    def _add_read_count(self, source,snv_idx, sample_idx, var=0, total=0):
    
    
        dat = self.file[f"read_counts/{source}"]
        dat['variant'][snv_idx,sample_idx] = var 
        dat['total'][snv_idx,sample_idx] = total
    
    def idx_by_label(self, label, index_name="SNV"):
        

        labs = self.file[f"{index_name}/label"]

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
                if snv_data == "log":
                    continue
                group_path = f"SNV/{snv_data}"
            
                self.file[group_path].resize((m,))
        else:
            m= self.file[f"SNV/label"].shape[0]

        if n:
            for sample_data in self.meta_tables:
                if sample_data == "log":
                    continue
                group_path = f"sample/{sample_data}"
                self.file[group_path].resize((n,))
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


    def _add_item(self, label:str, index_name, index_dict=None, cluster=None, data=None, overwrite=False):
        label = label.encode("utf-8")
        idx = self.idx_by_label(label, index_name, index_dict)

        if not idx:
            
            new_idx = self.file[f"{index_name}/label"].shape[0]
            if index_name == "SNV":
                self._resize_all(m =new_idx+1)
            else:
                self._resize_all(n=new_idx+1)
        else:
            if overwrite:
                new_idx = idx 
            else:
                print(f"{label} exists in {index_name} index_name, use overwrite=True to overwrite metadata.")
                return idx 
        
        self.file[f"{index_name}/label"][new_idx] = label
        self.file[f"{index_name}/index_map"][label] = new_idx
        if data:
            self.file[f"{index_name}/data"][new_idx] = data 
        
        if cluster:
            self.file[f"{index_name}/cluster"][new_idx] = cluster 


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




    def batch_add_snvs(self, labels, source_file=""):
      
        return self._batch_add_index(labels, index_name="SNV", source_file=source_file)
    


    def batch_add_samples(self, labels, source_file=""):
      
        return self._batch_add_index(labels, index_name="sample", source_file=source_file)
    

    def _batch_add_index(self, labels, index_name, source_file=""):
            index_dict = self._load_index(index_name)
            pre_size = len(index_dict)
      
            indices = []
            new = 0
            for lab in labels:
                if lab in index_dict:
                    indices.append(index_dict[lab])
                else:
                    next_idx = len(index_dict)
                    indices.append(next_idx)
                    index_dict[lab] = next_idx
                    new +=1
            indices = np.array(indices)
            labels  = np.array(labels)
            #sort the indices and labels concurrently for HDF5 fancy indexing
            labels, indices = labels[np.argsort(indices)].tolist(),  indices[np.argsort(indices)].tolist()
            if self.file[f"{index_name}/label"].shape[0 ]!= len(index_dict):
                if "SNV" == index_name:
    
                        self._resize_all(m=len(index_dict))
                
                if "sample" == index_name:
                    self._resize_all(n=len(index_dict))
        

            
            self.file[f"{index_name}/label"][indices] = labels
            self._save_index(index_dict, index_name)
            post_size = len(index_dict)
            if new > 0:

                self._update_index_log(index_name, new, pre_size, post_size, operation="add", source_file=source_file)

            if self.verbose:
                print(f"#{new} items added to {index_name} index")

            return index_dict
    

    def _update_index_log(self, index_name, num, pre_size, post_size, operation, source_file):
        source_file  = str(pathlib.Path(source_file).resolve())
        log =self.file[f"{index_name}/log"]
        current_size = log.shape[0]
        new_size = current_size + 1
        log.resize((new_size,))
        timestamp_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S").encode("utf-8")
        user = getpass.getuser()
        hostname = socket.gethostname().encode("utf8")
        log[current_size] = (timestamp_str, num, pre_size,post_size, operation,user, hostname, source_file)

    
    def get_snv_data(self, indices=None):

        return self._get_data(dataset_name="SNV/data", indices=indices)
    

    def get_snv_log(self):
        return self._get_data(dataset_name="SNV/log")
    
    def get_sample_log(self):
        return self._get_data(dataset_name="sample/log")
        
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
    

    def _update_value(self, indices, dataset_name, col, value):
            """
            Updates the value of a single column in a structured array
            indices: array-like object
            dataset_name: str dataset table name to update
            col: str: column name to update
            value: array-like or single-value to update the column of the specified indices,
                    if array, value must be sorted in the same order as indices
            """
            
            dataset = self.file[dataset_name]
            max_index = max(indices) + 1
            if dataset.shape[0] < max_index:
                self.file.close()
                raise IndexError("Invalid indices passed to update_data method, check indices and try again")

            # Read affected rows from HDF5 into memory
            temp_data = dataset[indices]  # Read only required rows
            
            # Update only the specified column
            temp_data[col] = value 

            # Write back only modified rows
            dataset[indices] = temp_data  

    def add_maf_files(self, fnames, **kwargs)  :
        """
        fnames: a list of paths to maf files
        """  
        for f in fnames:
            self.add_maf_file(f, **kwargs)


    def add_maf_file(self, fname, 
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

      


    
        snv_labels = maf["label"].tolist() 

        

        try: 

            
                snv_idx = self.batch_add_snvs(snv_labels, source_file=fname)

                #sort dataframe according to the newly assigned indices in DNAStream
                maf.loc[:,"snv_idx"] = maf["label"].map(snv_idx) 
                maf = maf.sort_values("snv_idx")
                indices = maf["snv_idx"].tolist()
                snv_data = maf[[val for _, val in column_dict.items()]]
              
    
                self.add_snv_data(indices, snv_data)


        except Exception as e:
       
            self.close()
            raise Exception(e)

    


    

    
    
    def close(self):
        self.file.close()

rc_pth = "/rsrch6/home/genetics/vanloolab/llweber/MPNST/scdna/read_counts"
ds = DNAStream(filename="temp.h5", verbose=True)
print(ds)
samples = [f"GEM2.2_PT_{i}" for i in range(2,6)] 
maf_files = [f"../data_summary/WGS_MUTATION/Somatic/2outof3_SNV/{sample}_SNVs_2outof3.maf" for sample in samples]
indices = ds.add_maf_files(maf_files)
print(ds.get_snv_log())
print(ds)
read_count_file = f"{rc_pth}/GEM2.2.csv"
ds.add_read_counts(read_count_file, source="scdna")
print(ds.get_snv_log())
print(ds.get_sample_log())
print(ds)
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

        






        