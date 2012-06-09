using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using utils;
using System.IO;

namespace ConsoleApplication1
{
    class Program
    {
        static void compute(string file, string fileo)
        {           
            double[] uk = utils.File.ReadData(file);
            int[] freq = { 17520,24*4};
            TimeSeries t = new TimeSeries(uk, freq);
            int[] freq1 = { 17520, 24 * 7*4 };
            TimeSeries t1 = new TimeSeries(uk, freq1);
            int[] freq2 = { 17520 };
            TimeSeries t2 = new TimeSeries(uk, freq2);
            Model m = new Model(t);

            m.Solve();
            m.Print();
            m.PrintShort();


            Model m1 = new Model(t1);

            m1.Solve();

            m1.PrintShort();


            Model m2 = new Model(t2);

            m2.Solve();

            m2.PrintShort();

            double[] d = new double[uk.Length];
            for (int i = 0; i < uk.Length; i++)
                d[i] = m.Eval(i);
            
            try
            {
                using (StreamWriter sw = new StreamWriter(fileo))
                {
                    for (int i = 0; i < uk.Length; i++)
                        sw.WriteLine(uk[i]+"\t"+d[i]);
                }
            }
            catch (Exception e)
            {
                // Let the user know what went wrong.
                Console.WriteLine("The file could not be written:");
                Console.WriteLine(e.Message);
            }

            
        }
        static void Main(string[] args)
        {
            Console.WriteLine("Global.use_values 0");
            Global.use_values = 0;
            //compute("c:/data/ukc.txt", "c:/data/a_c");
            Console.WriteLine("Global.use_values 0");
            Global.use_values = 1;
            compute("c:/data/ukc.txt", "c:/data/a_c");
          //  compute("c:/data/uk.txt", "c:/data/data");
             
        }
    }
}
