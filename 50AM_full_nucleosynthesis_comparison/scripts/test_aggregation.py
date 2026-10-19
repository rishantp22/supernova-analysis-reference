"""Scientific failure-mode reproduction on synthetic files, never real outputs."""
import subprocess,tempfile,unittest
from pathlib import Path
import pandas as pd

class AggregationMechanism(unittest.TestCase):
    def test_repeated_printed_time_merges_distinct_blocks(self):
        source=Path('/scratch/svierlin/prism-v02/sum_eoc_manifest_RP.sh').read_text()
        program=source[source.index('\nBEGIN {')+1:source.rindex("\n'")]
        with tempfile.TemporaryDirectory() as tmp:
            base=Path(tmp);(base/'7').mkdir()
            header='timestep {step} time(s) {time} temperature(GK) 1.0 density(g/cm^3) 1.0\n26 54 0.25\n'
            for times,expected in [([1000,2000],.25),([2000,2000],.5)]:
                (base/'7/xtime.txt').write_text('  Z A X(Z,A)\n'+''.join(header.format(step=k+1,time=t) for k,t in enumerate(times)))
                out=base/('case_'+str(times[0])+'.csv')
                subprocess.run(['awk','-v',f'BASE={base}','-v',f'OUTFILE={out}',program],input='7\n',text=True,check=True,capture_output=True)
                d=pd.read_csv(out,header=None,names=['id','category','step','time','Z','A','X'])
                value=d.loc[(d.category=='last')&(d.Z==26)&(d.A==54),'X'].item()
                self.assertEqual(value,expected)

if __name__=='__main__':unittest.main()
