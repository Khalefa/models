using System;
using System.Collections.Generic;
using System.Collections;
using System.Linq;
using System.Text;
using stllib;
using utils;

namespace ConsoleApplication1
{
    public class TimeSeries
    {
        public double []data;
        public int []freq;
        public TimeSeries() { }
        public TimeSeries(int n,int []freq)
        {
            data = new double[n];
            this.freq = freq;           
        }
        public TimeSeries(double []x, int []freq)
        {
            int n = x.Length;
            data = new double[n];
            for (int i = 0; i < n; i++) data[i] = x[i];
            this.freq = freq;
        }
        public int Length
        {
            get
            {
                return data.Length;
            }
        }

        public void Print()
        {
            if(freq==null)
                Console.WriteLine("Length" + data.Length + " Freq  null");
            else
            Console.WriteLine("Length" + data.Length + " Freq " + freq[0]);
        }

        public override string ToString()
        {
            string s="";
            for(int i=0;i<data.Length;i++)
                s=s+data[i]+ " ";
            return s;
        }
    }


}
