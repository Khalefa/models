using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ConsoleApplication1
{
    public class Range
    {
        public int s, e;
        public Range() { }
        public Range(int s, int e)
        { this.s = s; this.e = e; }

       public override string  ToString(){
           return "[" + s + "," + e + "]";
    }
    
    }
}
