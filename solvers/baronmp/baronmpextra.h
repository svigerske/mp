/*Adding 11 continuous and 0 integer variables with names.
- inf <= x <= inf
- inf <= y <= inf
1e-06 <= xx <= inf
0 <= yy <= 1
1e-06 <= c1 <= inf
0 <= c1_3_ <= 0
0 <= c2 <= 0
0 <= c5 <= 0
0 <= c5_3_ <= inf
0 <= c2_3_ <= 0
0 <= c5_5_ <= 0

c3: -inf <= x - 3 * y <= 3.3
c4 : -inf <= x - y <= 4
c5_4_ : c5_3_ = abs(y)

2 <= log(c1) <= inf

- inf <= ((x ^ 4)(y * x * 0.2)) <= 1.1
- inf <= -c5_3_ + ((x * x * -1)(yy * log(xx))) <= 0
c1 <= ((x * x))

z : maximize x - y
x - VISITOR 0.0.0 : not solved
9 continuous variables
2 integer variables
1 objective - linear
2 Linear constraints
1 Abs constraints
0 simplex iterations

------------WARNINGS------------*/