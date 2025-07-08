from gprofiler import GProfiler
import pandas as pd
df = pd.read_csv('upregulated_in_responders.csv')
df2 = pd.read_csv('downregulated.csv')
upregulated_genes = list(df['Description'])
downregulated_genes = list(df['Description'])
gp = GProfiler(return_dataframe=True)
results = gp.profile(organism="hsapiens", query=upregulated_genes)
results = gp.profile(organism ='hsapiens', query = downregulated_genes)
if results.empty:
    print("⚠️ No enrichment terms found. Some genes may be unrecognized.")
    print("Tip: try running with fewer, well-known gene symbols.")
else:
    results[['source', 'name', 'p_value', 'intersection_size']].to_csv('downregulated_analysis.csv')

#if results.empty:
 #   print("⚠️ No enrichment terms found. Some genes may be unrecognized.")
  #  print("Tip: try running with fewer, well-known gene symbols.")
#else:
 #   results[['source', 'name', 'p_value', 'intersection_size']].to_csv('upregulated_analysis_v2.csv')
