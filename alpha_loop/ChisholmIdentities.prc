#procedure ChisholmIdentities()
  label retry;
    repeat id d_(mu1?,mu2?)*d_(mu2?,mu3?) = d_(mu1,mu3);
    repeat id d_(mu1?,mu2?)*gamma(?a,mu1?,?b) = gamma(?a,mu2,?b);
    id once ifmatch->retry gamma(?a,mu?,mu?,?b)=gamma(?a,?b)*rat(4-2*ep,1);
    id once ifmatch->retry gamma(?a,mu?,mu1?,mu?,?b)=rat(-2+2*ep,1)*gamma(?a,mu1,?b);
    id once ifmatch->retry gamma(?a,mu?,mu1?,mu2?,mu?,?b)=
        +gamma(?a,?b)*4*d_(mu1,mu2)
        +rat(-2*ep,1)*gamma(?a,mu1,mu2,?b);
    id once ifmatch->retry gamma(?a,mu?,mu1?,mu2?,mu3?,mu?,?b)=
        -2*gamma(?a,mu3,mu2,mu1,?b)
        +rat(2*ep,1)*gamma(?a,mu1,mu2,mu3,?b);
    id once ifmatch->retry gamma(?a,mu?,mu1?,mu2?,mu3?,mu4?,mu?,?b) =
        +2*gamma(?a,mu4,mu1,mu2,mu3,?b)
        +2*gamma(?a,mu3,mu2,mu1,mu4,?b)
        +rat(-2*ep,1)*gamma(?a,mu1,mu2,mu3,mu4,?b);
    id once ifmatch->retry gamma(?a,mu?,mu1?,...,mu5?,mu?,?b) =
        +2*gamma(?a,mu5,mu1,...,mu4,?b)
        -2*gamma(?a,mu4,mu1,...,mu3,mu5,?b)
        -2*gamma(?a,mu3,mu2,mu1,mu4,mu5,?b)
        +rat(2*ep,1)*gamma(?a,mu1,...,mu5,?b);
    id once ifmatch->retry gamma(?a,mu?,mu1?,...,mu6?,mu?,?b) =
        +2*gamma(?a,mu6,mu1,...,mu5,?b)
        -2*gamma(?a,mu5,mu1,...,mu4,mu6,?b)
        +2*gamma(?a,mu4,mu1,...,mu3,mu5,mu6,?b)
        +2*gamma(?a,mu3,mu2,mu1,mu4,mu5,mu6,?b)
        +rat(-2*ep,1)*gamma(?a,mu1,...,mu6,?b);
    id once ifmatch->retry gamma(?a,mu?,mu1?,...,mu7?,mu?,?b) =
        +2*gamma(?a,mu7,mu1,...,mu6,?b)
        -2*gamma(?a,mu6,mu1,...,mu5,mu7,?b)
        +2*gamma(?a,mu5,mu1,...,mu4,mu6,mu7,?b)
        -2*gamma(?a,mu4,mu1,...,mu3,mu5,mu6,mu7,?b)
        -2*gamma(?a,mu3,mu2,mu1,mu4,mu5,mu6,mu7,?b)
        +rat(2*ep,1)*gamma(?a,mu1,...,mu7,?b);
* the gemueric fall-back case:
    id once ifmatch->retry gamma(?a,mu?,?b,mu1?,mu?,?c) =
        +2*gamma(?a,mu1,?b,?c) - gamma(?a,mu,?b,mu,mu1,?c);
#endprocedure