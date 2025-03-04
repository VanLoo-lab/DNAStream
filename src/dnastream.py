import pandas as pd 
import h5py 
"""
/
 ├── SNV/                     # Shared SNV index
 │   ├── labels               #short name chr:pos:ref:alt
     |-- data                 #dataframe structure containing quality scores, number of callers, etc
     |            #cluster ids
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
        self.h5file = h5py.File(filename, "a")  # Append mode (does not overwrite)

        # Create shared SNV index if not already present
        for snv_data in ["labels", "data", "cluster"]:
            group_path = [f"SNV/{snv_data}"]
            if group_path not in self.h5file:
                self.h5file.create_dataset(group_path, shape=(nsnvs,), maxshape=(None,), dtype=h5py.string_dtype())

        # Create group structure
        for modality in ["bulk", "lcm", "scdna"]:
                group_path = f"read_counts/{modality}"
                if group_path not in self.h5file:
                    self.h5file.create_group(group_path)
                for reads in ["variant", "total"]:
                    self.h5file.create_dataset(f"{group_path}/{reads}", shape=(nsnvs,), maxshape=(None,), dtype='i')

        
    


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

        






        