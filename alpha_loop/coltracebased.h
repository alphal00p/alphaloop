Symbol cOlNA,cOlNF;
AutoDeclare Index cOli=cOlNF, cOlii=cOlNF, cOlj=cOlNA, cOljj=cOlNA;
Set colF: cOli1,...,cOli80;
Set colA: cOlj1,...,cOlj80;
Set colFdum: cOlii1,...,cOlii80;
Set colAdum: cOljj1,...,cOljj80;
Tensor TCol,TrCol(cyclic),fCol(antisymmetric),ffCol;
CF countF, countA,color;
S aa,bb,cc,dd,ii,Tf,xx;


* write in terms of cf and ca 

#procedure tracebased
* small trick I learned from Ben. Allows for .sort and speeds up things.
    #define nloop "1"
******************
* Color evaluation
******************
    Multiply color(1);
    repeat id d_(cOlj1?,cOlj2?) * color(xx?) = color(xx*d_(cOlj1,cOlj2));
    repeat id d_(cOli1?,cOli2?) * color(xx?) = color(xx*d_(cOli1,cOli2));
    repeat id TCol(?aa)*color(xx?) = color(xx*TCol(?aa));
    repeat id fCol(cOlj1?,cOlj2?,cOlj3?)*color(xx?) = color(xx * fCol(cOlj1,cOlj2,cOlj3));
    repeat id ffCol(?aa)*color(xx?) = color(xx * ffCol(?aa));
    
    B+ color;
    .sort:color-prep;
    Keep Brackets;

* peform everything ones per color function    
    Argument color;
        Multiply countF(1);
        Multiply countA(1);
    EndArgument;
    B+ color;
    .sort:`nloop++';          
    Keep Brackets;
    #do loop=1,1
  
        Argument color;

* preparation
            id TCol(cOli1?,cOli2?) = d_(cOli1,cOli2);
            id TrCol = cOlNF;
            repeat id TCol(cOli1?,?aa,cOli2?)*TCol(cOli2?,?bb,cOli3?) = TCol(cOli1,?aa,?bb,cOli3);
            id  TCol(cOli1?,?aa,cOli1?) = TrCol(?aa);
            id  TrCol(cOlj1?) = 0;
            repeat id ffCol(cOlj1?,cOlj2?,cOlj3?,cOlj4?)*countA(ii?) = fCol(colAdum[ii], cOlj1,cOlj2)*fCol(colAdum[ii],cOlj3,cOlj4)*countA(ii+1);
            repeat id fCol(cOlj1?,cOlj2?,cOlj3?)= 2*i_*(TrCol(cOlj1,cOlj3,cOlj2)-TrCol(cOlj1,cOlj2,cOlj3));
* actual reduction
* 0. case: Outside of trace
            id ifmatch->end  TCol(cOli1?,?aa,cOlj1?,?bb,cOli2?)*TCol(cOli3?,?cc,cOlj1?,?dd,cOli4?) = Tf*( TCol(cOli1,?aa,?dd,cOli4)*TCol(cOli3,?cc,?bb,cOli2)- 1/cOlNF*TCol(cOli1,?aa,?bb,cOli2)*TCol(cOli3,?cc,?dd,cOli4));
            id ifmatch->end TCol(cOli1?,?aa,cOlj1?,?bb,cOlj1?,?cc,cOli2?) = Tf*( TCol(cOli1,?aa,?cc,cOli2)*TrCol(?bb)-1/cOlNF*TCol(cOli1,?aa,?bb,?cc,cOli2));
* 1. case: within trace
            id ifmatch->end TrCol(?aa,cOlj1?,?bb,cOlj1?,?cc)=Tf*(TrCol(?aa,?cc)*TrCol(?bb)-1/cOlNF*TrCol(?aa,?bb,?cc) );
            id ifmatch->end TrCol(?aa,cOlj1?,?bb)*TrCol(?cc,cOlj1?,?dd) =  Tf*( TrCol(?aa,?dd,?cc,?bb)- 1/cOlNF*TrCol(?aa,?bb)*TrCol(?cc,?dd));
* 2. case: outside and inside trace
            id ifmatch->end TCol(cOli1?,?aa,cOlj1?,?bb,cOli2?)*TrCol(?cc,cOlj1?,?dd) = Tf*( TCol(cOli1,?aa,?dd,?cc,?bb,cOli2) - 1/cOlNF*TCol(cOli1,?aa,?bb,cOli2)*TrCol(?cc,?dd));         
            goto end2;
            label end;    
            redefine loop "0";
            label end2;

            id TCol(cOli1?,cOli2?) = d_(cOli1,cOli2);
* TrCol() = TrCol(id) = cOlNF        
            id TrCol = cOlNF;  
        EndArgument;
            B+ color;
            .sort:reduction-`nloop++';
            Keep Brackets;
            
    #enddo
    Argument color;
        id countA(ii?) = 1;
        id countF(ii?) = 1;
    EndArgument;


#endprocedure