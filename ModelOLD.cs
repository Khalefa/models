using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using stllib;

namespace ConsoleApplication1
{
    class ModelOLD
    {
        //public enum Type { TREND, SEASON} ;
        public TimeSeries ts;
        public  ModelOLD trend;
        public ModelOLD season;
        public double[] values; // instead of trend, used when use_values
        
        public ModelOLD(TimeSeries ts)
        {
            this.ts = ts;
            trend = null;
            season = null;
            values = null;        
        }
        public ModelOLD()
        {
            ts = null;
            trend = null;
            season = null;
            values = null;        
        }

        private void decompose(int n, int freq,double[] trend_, double[] season_)
        {
            int swindow = 10 * n + 1;
            stllib.STL stl = new STL();
            double[] seaonal = new double[freq];
            int[] count = new int[freq];
            int limit= (int)Math.Ceiling((double)n / freq);
            if (Global.use_values == 1)
                values = new double[limit];
            else
                values = null;
            unsafe
            {
                double* y = (double*)utils.Memory.Alloc(sizeof(double) * n);
                double* t = (double*)utils.Memory.Alloc(sizeof(double) * n);
                double* s = (double*)utils.Memory.Alloc(sizeof(double) * n);
                for (int i = 0; i < n; i++) y[i] = ts.data[i];
                stl.compute(y, n, freq, swindow, t, s);
                for (int i = 0; i < n; i++)
                {
                    trend_[i] = t[i];
                }

                if (Global.use_values == 1)
                {
                    int k = 0;

                    for (int i = 0; i < limit; i++)
                    {
                        double t1 = 0;
                        double t2 = 0;
                        int l = 0;
                        for (int j = 0; j < freq; j++)
                        {
                            if (k == n) break;
                            l++;
                            k++;
                            t1 += ts.data[j + i * freq];
                            t2 += y[j + i * freq];
                        }
                        t1 /= l;
                        t2 /= l;
                        values[i] = (t2 + t1) / 2;
                    }
                }
                
                //make seaonality perfect
                for (int i = 0; i < n; i++)
                {
                    seaonal[i % freq] += s[i];
                    count[i % freq]++;
                }
                utils.Memory.Free(y);
                utils.Memory.Free(s);
                utils.Memory.Free(t);

            }
            for (int i = 0; i < freq; i++)
            {
                season_[i] = seaonal[i % freq] / count[i % freq];
            }           
            
            
        }
        public void Solve()
        {
            int n = ts.Length;
            if ((ts.freq == null) || (ts.freq[0] == 0))
            {
                if (Global.use_values == 0)
                {
                    double[] x = ChebyshevReg.Solve(ts.data);
                    values = new double[2];
                    values[0] = x[0];
                    values[1] = x[1];
                    season = null;
                    trend = null;
                }
                
            }
            else
            {
                int freq = ts.freq[0];

                double[] trend_ = new double[n];
                double[] season_ = new double[freq];

                decompose(n, freq, trend_, season_);        
                trend = new ModelOLD();
                trend.ts = new TimeSeries(trend_, null);
                season = new ModelOLD();
                int[] f;
                if (ts.freq.Length == 1) f = null;
                else
                f = new int[ts.freq.Length - 1];
                for (int i = 0; i < ts.freq.Length - 1; i++) f[i] = ts.freq[i + 1];
                season.ts = new TimeSeries(season_, f);
                trend.Solve();
                if(f!=null)
                    season.Solve();
            }
        }
        public double Eval(int x)
        {
            double t,s;
            t = s = 0;
            if (trend == null && season == null && values == null) return ts.data[x];
            
            if ((values != null) &&(Global.use_values==1) )
                t=values[x/ts.freq[0]];
            else
            if (trend!=null) t=trend.values[0]*x+trend.values[1];
            
            if(season!=null) s=season.Eval(x%ts.freq[0] );
            return s+t;
        }
        public double AvgError()
        {
            double error = 0;
            for (int i = 0; i < ts.Length; i++)
            {
                error += Math.Abs(ts.data[i] - Eval(i)) / ts.data[i];

            }

            return error/ts.Length*100;
        }
        public double MaxError()
        {
            double err = 0;
            double error = 0;
            for (int i = 0; i < ts.Length; i++)
            {
                err = Math.Abs(ts.data[i] - Eval(i)) / ts.data[i];
                if (err > error) error = err;

            }

            return error  * 100;
        }
        public int Size()
        {
            int size_ = 0;
            if (trend == null && season == null && values == null) size_ = ts.data.Length;
            else
            {
                if (values != null) size_ += values.Length;
                if (trend != null) size_ += trend.Size();
                if (season != null) size_ += season.Size();
            }
            return size_;
        }

        public void Print() {
            Console.WriteLine("{");
            ts.Print();
            Console.Write("Trend:");
            if (trend == null) Console.WriteLine("null");
            else trend.Print();
            Console.Write("Season:");
            if (season == null) Console.WriteLine("null");
            else season.Print();

            /*if (v != null)
            {
                Console.Write("v[" + v.Length + "]");
                Console.WriteLine(v[0] + "\t" + v[1]);
            }
            else Console.WriteLine("v:null");*/
            if (values != null)
            {
                Console.Write("values[" + values.Length + "]");
                Console.WriteLine(values[0]+ "\t"+ values[1]);
            }
            else Console.WriteLine("values:null");

            Console.WriteLine("Error:" + AvgError() + "% Size:" + Size());
            Console.WriteLine("}");
        }

        public void PrintShort()
        {

            Console.WriteLine("Error:" + AvgError() + "% " + "Max Error:" + MaxError() + "% Size:" + Size());
            
        }
    }

}
