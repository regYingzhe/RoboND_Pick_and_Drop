#!/usr/bin/env python

import numpy as np;
from numpy import array;
from sympy import symbols, cos, sin, pi, simplify, sqrt, atan2;
from sympy.matrices import Matrix;

### Create symols for joint variables
q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8');
d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8');
a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7');
alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7');


### KUKU KR210 ###
# DH Parameters
s = {
	alpha0:     0, a0:      0, d1:  0.75,
	alpha1: -pi/2, a1:   0.35, d2:     0, q2:q2 - pi/2,
	alpha2:     0, a2:   1.25, d3:     0, 
	alpha3: -pi/2, a3: -0.054, d4:   1.5, 
	alpha4:  pi/2, a4:      0, d5:     0, 
	alpha5: -pi/2, a5:      0, d6:     0, 
	alpha6:     0, a6:      0, d7: 0.303, q7:0 
}
print("d7 is: ", s[d7]);

##### Homogeneous Transforms
# bask_link to link1
T0_1 = Matrix([[            cos(q1),            -sin(q1),            0,                a0],
	           [sin(q1)*cos(alpha0), cos(q1)*cos(alpha0), -sin(alpha0), -sin(alpha0) * d1],
	           [sin(q1)*sin(alpha0), cos(q1)*sin(alpha0),  cos(alpha0),  cos(alpha0) * d1],
	           [                  0,                  0,             0,                 1]]);
T0_1 = T0_1.subs(s);

# link1 to link2
T1_2 = Matrix([[            cos(q2),            -sin(q2),            0,                a1],
			   [sin(q2)*cos(alpha1), cos(q2)*cos(alpha1), -sin(alpha1), -sin(alpha1) * d2],
			   [sin(q2)*sin(alpha1), cos(q2)*sin(alpha1),  cos(alpha1),  cos(alpha1) * d2],
			   [                  0,                  0,             0,                 1]]);
T1_2 = T1_2.subs(s);

# link2 to link3
T2_3 = Matrix([[            cos(q3),            -sin(q3),            0,                a2],
			   [sin(q3)*cos(alpha2), cos(q3)*cos(alpha2), -sin(alpha2), -sin(alpha2) * d3],
			   [sin(q3)*sin(alpha2), cos(q3)*sin(alpha2),  cos(alpha2),  cos(alpha2) * d3],
			   [                 0,                  0,             0,                 1]]);
T2_3 = T2_3.subs(s);

# link3 to link4
T3_4 = Matrix([[            cos(q4),            -sin(q4),            0,                a3],
	           [sin(q4)*cos(alpha3), cos(q4)*cos(alpha3), -sin(alpha3), -sin(alpha3) * d4],
	           [sin(q4)*sin(alpha3), cos(q4)*sin(alpha3),  cos(alpha3),  cos(alpha3) * d4],
	           [                 0,                  0,             0,                 1]]);
T3_4 = T3_4.subs(s);

# link4 to link5
T4_5 = Matrix([[            cos(q5),            -sin(q5),            0,                a4],
	           [sin(q5)*cos(alpha4), cos(q5)*cos(alpha4), -sin(alpha4), -sin(alpha4) * d5],
	           [sin(q5)*sin(alpha4), cos(q5)*sin(alpha4),  cos(alpha4),  cos(alpha4) * d5],
	           [                 0,                  0,             0,                 1]]);
T4_5 = T4_5.subs(s);

# link5 to link6
T5_6 = Matrix([[            cos(q6),            -sin(q6),            0,                a5],
	           [sin(q6)*cos(alpha5), cos(q6)*cos(alpha5), -sin(alpha5), -sin(alpha5) * d6],
	           [sin(q6)*sin(alpha5), cos(q6)*sin(alpha5),  cos(alpha5),  cos(alpha5) * d6],
	           [                 0,                  0,             0,                 1]]);
T5_6 = T5_6.subs(s);

# link6 to EE
T6_G = Matrix([[            cos(q7),            -sin(q7),            0,                a6],
	           [sin(q7)*cos(alpha6), cos(q7)*cos(alpha6), -sin(alpha6), -sin(alpha6) * d7],
	           [sin(q7)*sin(alpha6), cos(q7)*sin(alpha6),  cos(alpha6),  cos(alpha6) * d7],
	           [                 0,                  0,             0,                 1]]);
T6_G = T6_G.subs(s);

# Composition of homogeneous Transforms
T0_2 = simplify(T0_1 * T1_2);
T0_3 = simplify(T0_2 * T2_3);
T0_4 = simplify(T0_3 * T3_4);
T0_5 = simplify(T0_4 * T4_5);
T0_6 = simplify(T0_5 * T5_6);
T0_G = simplify(T0_6 * T6_G);

# correction Needed to Accout of Orientation Difference Between Definition of 
# Gripper_Link in URDF versus DH Convention

R_z = Matrix([[cos(np.pi),   -sin(np.pi),              0,     0],
	          [sin(np.pi),    cos(np.pi),              0,     0],
	          [         0,             0,              1,     0],
	          [         0,             0,              0,     1]]);
R_y = Matrix([[ cos(-np.pi/2),         0,  sin(-np.pi/2),     0],
	          [             0,         1,              0,     0],
	          [-sin(-np.pi/2),         0,  cos(-np.pi/2),     0],
	          [             0,         0,              0,     1]]);
R_corr = simplify(R_z * R_y);

T_total = simplify(T0_G * R_corr);


#### Numerically evaluate transforms (compare this with output of tf_echo!)
theta = {q1: 0.85, q2: 0.39, q3: -0.89, q4: -1.00, q5: -0.26, q6: 0.16};
# print("T0_1 = ", T0_1.evalf(subs=theta));
# print("T0_2 = ", T0_2.evalf(subs=theta));
# print("T0_3 = ", T0_3.evalf(subs=theta));
# print("T0_4 = ", T0_4.evalf(subs=theta));
# print("T0_5 = ", T0_5.evalf(subs=theta));
# print("T0_6 = ", T0_6.evalf(subs=theta));
# print("T0_G = ", T0_G.evelf(subs=theta));
print("T_total = ", T_total.evalf(subs=theta));
T_total = T_total.evalf(subs=theta);
extract_matrix = T_total[0:3, 0:3];
print("extract matrix is: ", extract_matrix);

### Find the Euler Angel, from base_link to EE
alpha = atan2(T_total[1,0], T_total[0,0])
beta = atan2(-T_total[2,0], sqrt(T_total[0,0] ** 2 + T_total[1,0] ** 2))
gamma = atan2(T_total[2,1], T_total[2,2])
print("alpha is ", alpha);
print("beta is", beta);
print("gamma is", gamma);
roll = gamma;
pitch = beta;
yaw = alpha;
R = symbols('R');
Y = symbols('Y');
P = symbols('P');

### Build Rotation Matrix extrinsic
RotYaw_z = Matrix([[cos(Y), -sin(Y), 0],
	               [sin(Y), cos(Y),  0],
	               [0,           0,  1]]);
RotPitch_y = Matrix([[cos(P),  0, sin(P)],
					 [0,       1,      0],
					 [-sin(P), 0, cos(P)]]);
RotRoll_x = Matrix([[1,      0,      0],
					[0, cos(R),-sin(R)],
					[0, sin(R), cos(R)]]);
# RotYaw_z = RotYaw_z.row_join(Matrix([[0],[0],[0]])).col_join(Matrix([[0,0,0,1]]));
# RotPitch_y = RotPitch_y.row_join(Matrix([[0],[0],[0]])).col_join(Matrix([[0,0,0,1]]));
# RotRoll_x = RotRoll_x.row_join(Matrix([[0],[0],[0]])).col_join(Matrix([[0,0,0,1]]));
R_corr_extracted = R_corr[0:3,0:3];
Rrpy = simplify(RotYaw_z * RotPitch_y * RotRoll_x * R_corr_extracted);
euler = {Y: yaw, P: pitch, R: roll};
Rrpy = Rrpy.evalf(subs=euler);
print("Rrpy is: ", Rrpy);
nx = Rrpy[0,2];
ny = Rrpy[1,2];
nz = Rrpy[2,2];
px = T_total[0,3];
py = T_total[1,3];
pz = T_total[2,3];
print("nx, ny, nz, px, py, pz: ", nx, ny, nz, px, py, pz);

### Calculate Wrist center Poistion
dG = 0.303
wx = px - dG * nx;
wy = py - dG * ny;
wz = pz - dG * nz;
print("wx, wy, wz", wx, wy, wz);
theta1 = atan2(wy, wx);
side_c = s[a2];
wcx = wz - s[d1];
x2 = sqrt(wx**2 + wy**2) - s[a1];
side_b = sqrt(x2**2 + wcx**2);
side_a = sqrt(s[a3]**2 + s[d4]**2);
print("a,b,c Length is: ", side_a, side_b, side_c)
