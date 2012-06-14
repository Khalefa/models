using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace utils
{
    class File
    {
        public static double[] ReadData(string file, int size=int.MaxValue)
        {
            ArrayList data = new ArrayList();
            try
            {
                // Create an instance of StreamReader to read from a file.
                // The using statement also closes the StreamReader.
                int i = 0;
                using (StreamReader sr = new StreamReader(file))
                {
                    String line;
                    // Read and display lines from the file until the end of
                    // the file is reached.
                    while ((line = sr.ReadLine()) != null)
                    {
                        i++;
                        data.Add(double.Parse( line));
                        if (i > size) break;
                    }
                }
            }
            catch (Exception e)
            {
                // Let the user know what went wrong.
                Console.WriteLine("The file could not be read:");
                Console.WriteLine(e.Message);
            }
            
            double[] d=new double[data.Count];
            for(int i=0;i<data.Count;i++) d[i]=(double)data[i];
            return d;
        }
        public static void WriteData(string file, double [] data)
        {   
            try
            {
                using (StreamWriter sw = new StreamWriter(file))
                {
                    for(int i=0;i<data.Length;i++)
                    sw.WriteLine(data[i]);
                }
            }
            catch (Exception e)
            {
                // Let the user know what went wrong.
                Console.WriteLine("The file could not be written:");
                Console.WriteLine(e.Message);
            }

            
        }
    
        public static void WriteData(string file, int freq,double[] data, double[] t,double[] s,double[] e)
        {
            try
            {
                using (StreamWriter sw = new StreamWriter(file))
                {
                    for (int i = 0; i < data.Length; i++)
                        sw.WriteLine(data[i]+"\t"+t[i]+"\t"+s[i%freq]+"\t"+e[i]);
                }
            }
            catch (Exception ee)
            {
                // Let the user know what went wrong.
                Console.WriteLine("The file could not be written:");
                Console.WriteLine(ee.Message);
            }


        }
    }
}
