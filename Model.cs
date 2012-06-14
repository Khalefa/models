using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using stllib;
namespace ConsoleApplication1
{
    public enum ModelType
    {
        Explicit, Implicit, Trend
    };
    class ModelBase
    {
    }
    class Model : ModelBase
    {
        protected TimeSeries ts;
        protected ModelType type;

        protected double[] errors;
        double[] values;
        double[] trend;
        protected Model seasonal;
        int limit;
        int freq;
        private void decompose(int n, int freq, double[] season_)
        {
            int swindow = 10 * n + 1;
            stllib.STL stl = new STL();
            double[] seaonal = new double[freq];
            int[] count = new int[freq];
            int limit = (int)Math.Ceiling((double)n / freq);
            this.limit = limit;
            this.freq = freq;
            
            unsafe
            {
                double* y = (double*)utils.Memory.Alloc(sizeof(double) * n);
                double* t = (double*)utils.Memory.Alloc(sizeof(double) * n);
                double* s = (double*)utils.Memory.Alloc(sizeof(double) * n);
                for (int i = 0; i < n; i++) y[i] = ts.data[i];
                stl.compute(y, n, freq, swindow, t, s);
                int k = 0;
                if (type == ModelType.Explicit || type == ModelType.Implicit)
                {
                    values = new double[limit];
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
                if (type == ModelType.Implicit)
                {
                    double[] x = ChebyshevReg.Solve(values);
                    values = new double[2];
                    values[0] = x[0];
                    values[1] = x[1];
                }
                else if (type == ModelType.Trend)
                {
                    double[] t_ = new double[n];
                    for (int i = 0; i < n; i++)
                    {
                        t_[i] = t[i];
                    }

                    double[] x = ChebyshevReg.Solve(t_);
                    values = new double[2];
                    values[0] = x[0];
                    values[1] = x[1];
                }
                //make seaonality perfect
                for (int i = 0; i < n; i++)
                {
                    seaonal[i % freq] += s[i];
                    count[i % freq]++;
                }
                for (int i = 0; i < freq; i++)
                {
                    season_[i] = seaonal[i % freq] / count[i % freq];
                }

                trend = new double[n];
                errors = new double[n];
                for (int i = 0; i < n; i++)
                {
                    trend[i] = t[i];
                    errors[i] = ts.data[i] - (t[i] + s[i % freq]);
                }

                utils.Memory.Free(y);
                utils.Memory.Free(s);
                utils.Memory.Free(t);

            }

        }
        public Model(TimeSeries ts)
        {
            this.ts = ts;
            seasonal = null;
            values = null;
            errors = null;

            ModelType best = ModelType.Trend;
            double min = double.MaxValue;
            foreach (ModelType t in Enum.GetValues(typeof(ModelType)))
            {
                this.type = t;                
                this.Solve();
                this.CalcError();
                double rt = this.Size() * Error(0.9);
                if (min > rt) { min = rt; best = t; }
            }

            this.type = best;
            this.Solve();
            this.CalcError();
        }

        public int countError(double error)
        {
            int count = 0;
            if (errors == null) ; CalcError();
            
            for (int i = 0; i < ts.Length; i++)
            {
               if( errors[i]<error) count++;
            }

            return count;
        }
        public Model()
        {
            ts = null;
            seasonal = null;
            values = null;
            errors = null;
        }
        public void Solve()
        {
            int n = ts.Length;
            if ((ts.freq == null))
            {
                // ts data is used
                int x = 9;
                x++;
            }
            else
            {
                int freq = ts.freq[0];
                int l = 0;
                while (freq > ts.Length) { freq = ts.freq[l++]; }
                if (freq == 0)
                {
                    //use regression
                    values = ChebyshevReg.Solve(ts.data);
                    seasonal = null;
                }
                else
                {
                    double[] season_ = new double[freq];

                    decompose(n, freq, season_);
                    this.freq = freq;
                    seasonal = new Model();
                    int[] f;
                    if (ts.freq.Length == 1) f = null;
                    else
                    {
                        f = new int[ts.freq.Length - 1];
                        for (int i = 0; i < ts.freq.Length - 1; i++) f[i] = ts.freq[i + 1];
                    }
                    seasonal.ts = new TimeSeries(season_, f);

                    if (f != null)
                        seasonal.Solve();
                }
            }
        }


        public double Eval(int x)
        {
            if (seasonal == null && values == null) return ts.data[x];
            if (seasonal == null) return values[0] * x + values[1];
            double t, s;
            t = s = 0;
            if (type == ModelType.Explicit)
                t = values[x / freq];
            else if (type == ModelType.Implicit)
                t = values[0] * x / freq + values[1];
            else t = values[0] * x + values[1];

            s = seasonal.Eval(x % freq);
            return s + t;
        }
        private void CalcError()
        {
            errors = new double[ts.Length];
            for (int i = 0; i < ts.Length; i++)
            {
                errors[i] = Math.Abs(ts.data[i] - Eval(i)) / ts.data[i] * 100;
            }
            Array.Sort(errors);
        }
        public double AvgError()
        {
            if (errors == null) ; CalcError();
            double error = 0;
            for (int i = 0; i < ts.Length; i++)
            {
                error += errors[i];
            }

            return error / ts.Length;
        }
        public double MaxError()
        {
            if (errors == null) ; CalcError();

            return errors[ts.Length - 1];
        }
        public double Error(double x)
        {
            if (errors == null) ; CalcError();

            return errors[(int)(x * (ts.Length - 1))];
        }
        public virtual int Size()
        {
            if (seasonal == null && values == null) return ts.data.Length;
            if (seasonal == null) return 2;
            return values.Length + seasonal.Size();
        }
        public void Print()
        {
            Console.WriteLine("{");
            ts.Print();
            Console.Write("Season:");
            if (seasonal == null) Console.WriteLine("null");
            else seasonal.Print();
            if (values != null)
            {
                Console.Write("values[" + values.Length + "]");
                Console.WriteLine(values[0] + "\t" + values[1]);
            }
            else Console.WriteLine("values:null");

            Console.WriteLine("Error:" + AvgError() + "% Size:" + Size());
            Console.WriteLine("}");
        }
        public void PrintShort()
        {
            Console.WriteLine(ts.Length + "\t" + String.Format("{0:00.00}", Error(0.5)) + "% " + "\t" + String.Format("{0:00.00}", Error(0.9)) + "% " + "\t" + String.Format("{0:00.00}", MaxError()) + "%\t" + Size());
        }
        public void Save()
        {
            utils.File.WriteData("c:/data/d1", ts.freq[0], ts.data, trend, seasonal.ts.data, errors);
        }
    }

    class ModelSet
    {
        TimeSeries ts;
        int pieces;
        ModelType type;
        Model[] models;
        double[] errors;
        public ModelSet(TimeSeries ts, ModelType type, int pieces)
        {
            this.pieces = pieces;
            this.ts = ts;
            this.type = type;
            this.models = new Model[pieces];

            for (int i = 0; i < pieces; i++)
            {
                double[] u = new double[ts.Length / pieces];
                for (int j = 0; j < ts.Length / pieces; j++) u[j] = ts.data[i * ts.Length / pieces + j];
                TimeSeries t = new TimeSeries(u, ts.freq);
                Model m = new Model(t);
                models[i] = m;
            }

        }
        public void Solve()
        {
            for (int i = 0; i < pieces; i++)
            {
                models[i].Solve();
            }
        }
        public double Eval(int x)
        {
            int n = ts.Length / pieces;
            int l = x / n;
            return models[l].Eval(x % n);
        }
        public int Size()
        {
            int size = 0;
            for (int i = 0; i < pieces; i++)
            {
                size += models[i].Size();
            }
            return size;
        }
        private void CalcError()
        {
            errors = new double[ts.Length];
            for (int i = 0; i < ts.Length; i++)
            {
                errors[i] = Math.Abs(ts.data[i] - Eval(i)) / ts.data[i] * 100;
            }
            Array.Sort(errors);
        }
        public double AvgError()
        {
            if (errors == null) ; CalcError();
            double error = 0;
            for (int i = 0; i < ts.Length; i++)
            {
                error += errors[i];
            }

            return error / ts.Length;
        }
        public double MaxError()
        {
            if (errors == null) ; CalcError();

            return errors[ts.Length - 1];
        }
        public double Error(double x)
        {
            if (errors == null) ; CalcError();

            return errors[(int)(x * (ts.Length - 1))];
        }
        public void PrintShort()
        {
            Console.WriteLine(ts.Length + "\t" + String.Format("{0:00.00}", Error(0.5)) + "% " + "\t" + String.Format("{0:00.00}", Error(0.9)) + "% " + "\t" + String.Format("{0:00.00}", MaxError()) + "%\t" + Size());
        }
    }

    class ModelTree : Model
    {
        ModelTree[] children;
        Range range;
        double error;
        double[] errors;
        #region old
        /*  ModelTree[] findChildren(int len, ref double max_error)
        {
            max_error = double.MinValue;
            int kk = (int)Math.Floor((double)ts.Length / len);

            ModelTree[] c = new ModelTree[kk];
            int x = 0;

            for (int i = 0; i < kk - 1; i++)
            {
                if (i == kk - 2) len = ts.Length - x;

                double[] u = new double[len];
                for (int j = 0; j < len; j++) u[j] = ts.data[x++];
                TimeSeries t = new TimeSeries(u, ts.freq);
                ModelTree m = new ModelTree(t,  len * i);
                if (m.error > max_error) max_error = m.error;
                if (max_error > 0.9 * error) return null;
                c[i] = m;
            }

            if (max_error > 0.9 * error) return null;
            return c;
        }
        void setChildren( int len)
        {
            int kk = (int)Math.Floor((double)ts.Length / len);

            children = new ModelTree[kk];
            int x = 0;
            int l = ts.Length / kk;
            for (int i = 0; i < kk; i++)
            {
                if (i == kk - 1) l = ts.Length - x;
                double[] u = new double[l];
                for (int j = 0; j < l; j++) u[j] = ts.data[x++];
                TimeSeries t = new TimeSeries(u, ts.freq);
                ModelTree m = new ModelTree(t, l * i);
                children[i] = m;
            }
        }
        */

        /*  //divide the time series into k
            /*range = new Range(start + 0, start + ts.Length - 1);
            error = Error(0.9);
            if ((ts.Length < 100)) children = null;
            else
            {
                children = null;
                double max_error = 0;
                double min_max_error = double.MaxValue;
                int best_k = ts.freq[0];

                int i = 0;
                for (i = 0; i < ts.freq.Length; i++)
                {
                    int f = ts.freq[i];
                    if (f == 0) continue;
                    //       for (; ; )
                    {
                        ModelTree[] c = findChildren( f * 4, ref max_error);
                        // if (max_error < min_max_error) {
                        best_k = f * 8 / 2;
                        /*if (c != null)
                        {
                            children = new ModelTree[c.Length];
                            for (int m = 0; m < c.Length; m++)
                                children[m] = c[m];
                            break;
                        }*/
        /*      }
          }
          if (children == null) setChildren( best_k);

      }*/
        #endregion
        Model[] getModels(int len, double error, ref int length)
        {
            length = 0;
            int kk = (int)Math.Floor((double)ts.Length / len);

            Model[] c = new Model[kk];
            int x = 0;

            for (int i = 0; i < kk; i++)
            {
                if (i == kk - 1) len = ts.Length - x;

                double[] u = new double[len];
                for (int j = 0; j < len; j++) u[j] = ts.data[x++];
                TimeSeries t = new TimeSeries(u, ts.freq);
                Model m = new Model(t);
                length += m.countError(error);
                c[i] = m;
            }
            return c;
        }

        void setModels(int len)
        {
            
            int kk = (int)Math.Floor((double)ts.Length / len);

            children = new ModelTree[kk];
            int x = 0;

            for (int i = 0; i < kk; i++)
            {
                if (i == kk - 1) len = ts.Length - x;
                int s = x;
                double[] u = new double[len];
                for (int j = 0; j < len; j++) u[j] = ts.data[x++];
                TimeSeries t = new TimeSeries(u, ts.freq);
                ModelTree m = new ModelTree(t,null,s);
                children[i] = m;
            }
            
        }
        
        public void BuildTree()
        {
            //start with the lowest frequncy greater than zero
            int i = 0;
            for (int l = errors.Length - 1; l >= 0; l--)
            {
                double current_error = errors[l];
                int length = ts.Length;
                int len = 0;
                int best_f = 0;
                for (i = 0; i < ts.freq.Length; i++)
                {
                    int f = ts.freq[i];
                    if (f == 0) continue;
                    for (; ; )
                    {
                        Model[] c = getModels(f, current_error, ref len);
                        if (len >= 0.99 * length)
                        {
                            best_f = f;
                        }
                        else break;
                        f = f * 2;
                        if (len == 0) break;
                        if (c.Length == 1) break;
                        if ((i + 1 < ts.freq.Length) && (f > ts.freq[i + 1])) break;
                    }
                    Console.WriteLine(best_f);
                    if (best_f != 0) setModels(best_f);
                }
            }
           
        }

        public ModelTree(TimeSeries ts, double[] errors = null, int start = 0)
            : base(ts)
        {
            this.range = new Range(start, start + ts.Length - 1);
            this.errors = errors;
            this.children = null;
        }
        public override int Size()
        {
            int size = base.Size();
            size += 2;

            return size;
        }
        public string ToString(int indent = 0)
        {
            string h = indent + " ";
            for (int i = 0; i < indent; i++) h += "\t";
            int k = 0;
            if (children != null)
            {
                k = children.Length;
            }
            string s = h + "Model Error:" + Error(0.9) + "%" + range.ToString() + " " + Size() + " " + type + " " + k + "\n";
            if (children == null)
            {
                return s;
            }

            for (int i = 0; i < k; i++)
            {
                s += children[i].ToString(indent + 1);
            }
            return s;
        }
    }

}
