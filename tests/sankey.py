"""Aggregate Sankey contracts: retain unknowns and reject incomplete partitions."""
import copy,importlib.util,unittest
from pathlib import Path
spec=importlib.util.spec_from_file_location('sankey_plot',Path(__file__).resolve().parents[1]/'diagnostics/sankey/plot.py')
sankey=importlib.util.module_from_spec(spec);spec.loader.exec_module(sankey)
class SankeyCounts(unittest.TestCase):
 def setUp(self):
  self.rows=[dict(category=s,group='P1',cohort='primary_phenotypes',variable='first_hour_support',n=n,denominator_n='3',cell_status='observed',partition_complete='TRUE') for s,n in [('UNKNOWN','2'),('IMV','1')]]
  self.den=[dict(cohort='primary_phenotypes',group='P1',n='3'),dict(cohort='primary_phenotypes',group='I',n='0')]
  self.avail=[dict(cohort='primary_phenotypes',group=g,variable='first_hour_support',status=v) for g,v in [('P1','calculated'),('I','structural_zero_group')]]
 def call(self):return sankey.extract_counts(self.rows,self.den,self.avail,'phenotype')
 def test_unknown_and_zero_group(self):
  flows,den=self.call();self.assertEqual(flows,{('UNKNOWN','P1'):2,('IMV','P1'):1});self.assertEqual(den['I'],0)
 def test_no_overall_double_count(self):
  r=copy.deepcopy(self.rows[0]);r['group']='Overall';r['n']='100';self.rows.append(r)
  self.assertEqual(sum(self.call()[0].values()),3)
 def test_missing_counts_rejected(self):
  self.rows.pop()
  with self.assertRaises(ValueError):self.call()
 def test_suppression_rejected(self):
  self.rows[0]['cell_status']='suppressed_partition'
  with self.assertRaises(ValueError):self.call()
 def test_duplicate_cells_rejected(self):
  self.rows.append(copy.deepcopy(self.rows[0]))
  with self.assertRaises(ValueError):self.call()
 def test_unavailable_not_zero(self):
  self.avail[0]['status']='not_supplied_or_not_derived'
  with self.assertRaises(ValueError):self.call()
 def test_local_labels_kept_local(self):
  for table in [self.rows,self.den,self.avail]:
   for r in table:r['cohort']='provisional_local_clusters';r['group']='1' if r['group']=='P1' else '2'
  flows,den=sankey.extract_counts(self.rows,self.den,self.avail,'local_cluster')
  self.assertEqual(flows[('UNKNOWN','1')],2);self.assertEqual(den['2'],0)
if __name__=='__main__':unittest.main()
