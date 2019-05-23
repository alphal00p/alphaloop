import vectors; import math; sqrt = math.sqrt
euclidean_2_to_2_points = [
    vectors.LorentzVector([3.,0.,0.,sqrt(14.)])*(1./sqrt(5.)),
    vectors.LorentzVector([2.,0.,0.,-sqrt(14.)])*(1./sqrt(5.)),
    vectors.LorentzVector([-3.,sqrt(163./70.),sqrt(489./40.),-23.*sqrt(1./56.)])*(1./sqrt(5.)),
    vectors.LorentzVector([-2.,-sqrt(163/70.),-sqrt(489./40.),23.*sqrt(1./56.)])*(1./sqrt(5.))
]

physical_2_to_2_points = [
    vectors.LorentzVector([29.,0.,0.,sqrt(721.)])*(1./sqrt(120.)),
    vectors.LorentzVector([31.,0.,0.,-sqrt(721.)])*(1./sqrt(120.)),
    vectors.LorentzVector([-29.,4.*sqrt(366./103.),6.*sqrt(854./103.),-43.*sqrt(7./103.)])*(1./sqrt(120.)),
    vectors.LorentzVector([-31.,-4.*sqrt(366./103.),-6.*sqrt(854./103.),43.*sqrt(7./103.)])*(1./sqrt(120.))
]


def test_PS_point(PS):
    print("M1=",PS[0].square())
    print("M2=",PS[1].square())
    print("M3=",PS[2].square())
    print("M4=",PS[3].square())
    print("s=",(PS[0]+PS[1]).square())
    print("t=",(PS[0]+PS[2]).square())
    print("Sum=",(PS[0]+PS[1]+PS[2]+PS[3]))

test_PS_point(euclidean_2_to_2_points)
test_PS_point(physical_2_to_2_points)
