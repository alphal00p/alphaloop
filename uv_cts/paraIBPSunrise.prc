#procedure paraIBPSunrise() 
    #define nloop "1"
    #do loop=1,1
    #$doloop = 0;
* zero sectors
        id topo2(n1?{<1},n2?{<1},n3?{<1}) = 0;

        id topo2(n1?{<1},n2?{<1},n3?) = 0;

        id topo2(n1?{<1},n2?,n3?{<1}) = 0;

        id topo2(n1?,n2?{<1},n3?{<1}) = 0;

* sector mappings   

        id topo2(n1?{>0}, n2?{<1}, n3?{>0}) =
           topo2(n2, n3, n1);

        id topo2(n1?{>0}, n2?{>0}, n3?{<1}) =
            topo2(n3, n2, n1);

* actual reduction
     
        id ifmatch->end topo2(n1?{<0}, n2?{>0}, n3?{>1}) = 
            topo2(1 + n1, -1 + n2, n3)*rat(1,1) + 
            topo2(1 + n1, n2, -1 + n3)*rat((-3 + d - 3*n1),2*(-1 + n3)) + 
            topo2(2 + n1, -1 + n2, -1 + n3)*rat((-1 - n1),2*(-1 + n3)) + 
            topo2(2 + n1, n2, -2 + n3)*rat((1 + n1),2*(-1 + n3)) - 
            topo2(2 + n1, n2, -1 + n3)*rat(3*mUV^2*(1 + n1),2*(-1 + n3));

        id ifmatch->end topo2(n1?{<0}, n2?{>1}, n3?{>0}) = 
            topo2(1 + n1, -1 + n2, n3)*rat((-3 + d - 3*n1),2*(-1 + n2)) + 
             topo2(1 + n1, n2, -1 + n3)*rat(1,1) + 
            topo2(2 + n1, -2 + n2, n3)*rat((1 + n1),2*(-1 + n2)) + 
            topo2(2 + n1, -1 + n2, -1 + n3)*rat((-1 - n1),2*(-1 + n2)) - 
            topo2(2 + n1, -1 + n2, n3)*rat(3*mUV^2*(1 + n1),2*(-1 + n2));

        id ifmatch->end topo2(n1?{<1}, n2?{>0}, n3?{>1})  = 
            topo2(n1, n2, -1 + n3)*rat((2 + d - n1 - 2*n3),2*mUV^2*(-1 + n3)) + 
            topo2(1 + n1, -1 + n2, -1 + n3)*rat(n1,2*mUV^2*(-1 + n3)) - 
            topo2(1 + n1, n2, -2 + n3)*rat(n1,2*mUV^2*(-1 + n3)) - 
            topo2(1 + n1, n2, -1 + n3)*rat(n1,2*(-1 + n3));

        id ifmatch->end topo2(n1?{<1}, n2?{>1}, n3?{>0})  = 
            topo2(n1, -1 + n2, n3)*rat((2 + d - n1 - 2*n2),2*mUV^2*(-1 + n2)) - 
            topo2(1 + n1, -2 + n2, n3)*rat(n1,2*mUV^2*(-1 + n2)) + 
            topo2(1 + n1, -1 + n2, -1 + n3)*rat(n1,2*mUV^2*(-1 + n2)) - 
            topo2(1 + n1, -1 + n2, n3)*rat(n1,2*(-1 + n2));

        id ifmatch->end topo2(n1?{<0}, 1, 1)  = 
            topo2(1 + n1, 0, 1)*rat(n1,-2 + d - n1) - 
            topo2(1 + n1, 1, 0)*rat(n1,-2 + d - n1) + 
            topo2(1 + n1, 1, 1)*rat(mUV^2*(-3 + d - 2*n1),-2 + d - n1) + 
            topo2(1 + n1, 2, 0)*rat(2*mUV^2,-2 + d - n1) - 
            topo2(2 + n1, 0, 1)*rat(mUV^2*(1 + n1),2 - d + n1) - 
            topo2(2 + n1, 1, 0)*rat(mUV^2*(1 + n1),-2 + d - n1) - 
            topo2(2 + n1, 1, 1)*rat(3*mUV^4*(1 + n1),-2 + d - n1);

        id ifmatch->end topo2(n1?{>0}, n2?{>0}, n3?{>1})  = 
            topo2(-1 + n1, 1 + n2, -1 + n3)*rat(n2,3*mUV^2*(-1 + n3)) + 
            topo2(n1, n2, -1 + n3)*rat((3 + d - 3*n3),3*mUV^2*(-1 + n3)) - 
            topo2(n1, 1 + n2, -2 + n3)*rat(n2,3*mUV^2*(-1 + n3)) + 
            topo2(1 + n1, -1 + n2, -1 + n3)*rat(n1,3*mUV^2*(-1 + n3)) - 
            topo2(1 + n1, n2, -2 + n3)*rat(n1,3*mUV^2*(-1 + n3));          

        id ifmatch->end topo2(n1?{>0}, n2?{>1}, n3?{>0}) = 
            topo2(-1 + n1, n2, n3)*rat(1,3*mUV^2) + 
            topo2(n1, -1 + n2, n3)*rat((3 + d - 3*n2),3*mUV^2*(-1 + n2)) - 
            topo2(n1, n2, -1 + n3)*rat(1,3*mUV^2) - 
            topo2(1 + n1, -2 + n2, n3)*rat(2*n1,3*mUV^2*(-1 + n2)) + 
            topo2(1 + n1, -1 + n2, -1 + n3)*rat(2*n1,3*mUV^2*(-1 + n2)); 
        
        id ifmatch->end topo2(n1?{>1}, n2?{>0}, n3?{>0}) = 
            topo2(-2 + n1, 1 + n2, n3)*rat(-2*n2,3*mUV^2*(-1 + n1)) + 
            topo2(-1 + n1, n2, n3)*rat((3 + d - 3*n1),3*mUV^2*(-1 + n1)) + 
            topo2(-1 + n1, 1 + n2, -1 + n3)*rat(2*n2,3*mUV^2*(-1 + n1)) + 
            topo2(n1, -1 + n2, n3)*rat(1,3*mUV^2) - 
            topo2(n1, n2, -1 + n3)*rat(1,3*mUV^2);

        id ifmatch->end topo2(n1?{>0}, n2?{>0}, n3?{>1})  = 
            topo2(n1, 1 + n2, -1 + n3)*rat(n2,-1 + n3) + 
            topo2(1 + n1, -1 + n2, n3)*rat(1,1) - 
            topo2(1 + n1, 1 + n2, -2 + n3)*rat(n2,-1 + n3) + 
            topo2(2 + n1, -1 + n2, -1 + n3)*rat((-1 - n1),-1 + n3) + 
            topo2(2 + n1, n2, -2 + n3)*rat((1 + n1),-1 + n3); 

        label end;

        id topo2(0, 1, 1) =tp2P011;
        id topo2(1, 1, 1) = tp2P111;

        if (count(topo2, 1));
            $doloop = 1;
*            Print "%t";
        endif;
    
        ModuleOption maximum, $doloop;
        .sort:reduction-`nloop++';
        #if `$doloop'
            #redefine loop "0"
        #endif
    #enddo
* I want the MIs to be topo2(0,2,2) and topo2(2,2,1)
    id tp2P011=mi1L2*rat(4*mUV^4,d^2 - 4*d + 4);
    id tp2P111= rat(9*mUV^4,d^2 - 5*d + 6)*(mi2L2-mi1L2*rat(-1,3*mUV^2) );
#endprocedure