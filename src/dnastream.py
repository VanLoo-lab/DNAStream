import pandas as pd 
import h5py 
import numpy as np
"""
/
 ├── SNV/                     # Shared SNV index
 │   ├── labels               #short name chr:pos:ref:alt
     |-- data                 #dataframe structure containing quality scores, number of callers, etc
     |-- cluster
 ├── sample/                     # Shared SNV index
 │   ├── labels               #short name chr:pos:ref:alt
     |-- data                 #dataframe structure containing bam file path, sample code
     |-- cluster
 ├── read_counts/               # Read count matrices
 │   ├── bulk/                  # Bulk sequencing read counts
 │   │   ├── variant_reads       # SNVs x Samples (variant read counts)
 │   │   ├── total_reads         # SNVs x Samples (total read counts)
 │   ├── lcm/                    # LCM sequencing read counts
 │   │   ├── variant_reads       
 │   │   ├── total_reads         
 │   ├── scDNA/                  # scDNA-seq read counts
 │   │   ├── variant_reads       
 │   │   ├── total_reads         
 ├── metadata/                   # Metadata storage
 │   ├── sample_info              # Sample IDs
 │   ├── processing_parameters
"""



class DNAStream:
      def __init__(self, filename, nsnvs=0):
        """Initialize HDF5 storage."""
        self.filename = filename
        self.file = h5py.File(filename, "a")  # Append mode (does not overwrite)

        # snv_grp = self.self.h5file.create_group("SNV")
        # snv_grp = snv_grp.attrs['label'] = np.array(shape=0, dtype=str)
        # snv_grp = snv_grp.attrs['cluster']= np.array(shape=0, dtype=np.int32)
        self.meta_tables = ["data", "label", "cluster"]
        self.modalities = ["bulk", "lcm", "scdna"]
        self.indices = ["SNV", "sample"]

        # Create shared SNV index if not already present
        for index in self.indices:
            for data in self.meta_tables:
                group_path = [f"{index}/{data}"]
                if group_path not in self.file:
                    self.file.create_dataset(group_path, shape=(0,), maxshape=(None,), dtype=h5py.string_dtype())
          

        # Create group structure
        for modality in self.modalities:
                group_path = f"read_counts/{modality}"
                if group_path not in self.file:
                    self.file.create_group(group_path)
                for reads in ["variant", "total"]:
                    self.file.create_dataset(f"{group_path}/{reads}", shape=(0,0), maxshape=(None,), dtype='i')
                    self.file[f"{group_path}/{reads}"].dims[0].label = "SNV"
                    self.file[f"{group_path}/{reads}"].dims[1].label = "sample"

        
    

        def add_read_counts(self, fname, source):
            rc = pd.read_csv(fname)
            for _, row in rc.iterrows():
                self._add_read_count(source, row[0], row[1], row[3])

           

        
        def idx_by_label(self, label, index="SNV"):
            indices = np.where(self.file[f"{index}/label"]==label)[0]
            if len(indices) ==0:
                return None 
            elif len(indices) == 1:
                return indices[0]
            else:
                raise ValueError("SNV label is associated with multiple indices and could not be added.")
        

        def label_by_idx(self, idx, index="SNV"):
   
            return self.file[f"{index}/label"][idx]


        def _resize_all(self, m=None,n=None ):

            if m:
                for snv_data in self.meta_tables:
                    group_path = [f"SNV/{snv_data}"]
             
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
            return self.add_index(label, index="SNV", cluster=cluster, data=data, overwrite=overwrite)
        

        def add_sample(self, label, cluster=None, data=None, overwrite=False):
            return self.add_index(label, index="sample", cluster=cluster, data=data, overwrite=overwrite)


        def add_index(self, label, index, cluster=None, data=None, overwrite=False):
            idx = self.get_idx_by_label(self, label, index)

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
            
            if data:
                self.file[f"{index}/data"][new_idx] = data 
            
            if cluster:
                self.file[f"{index}/cluster"][new_idx] = cluster 
            else:
                self.file[f"{index}/cluster"][new_idx] = np.nan

            return new_idx





        def __add_read_count(self, source, snv_label, sample_label, var=0, total=0):
            snv_idx = self.add_snv(snv_label)
            sample_idx = self.add_sample(sample_label)
        
            dat = self.file[f"read_counts/{source}"]
            dat['variant'][snv_idx,sample_idx] = var 
            dat['total'][snv_idx,sample_idx] = total


                  
                  
                  
             

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

        






        