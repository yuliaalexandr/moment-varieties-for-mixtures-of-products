R = QQ[x12,x13,x14,x15,x16,x21,x23,x24,x25,x26,x31,x32,x34,x35,x36,x41,x42,x43,x45,x46,x51,x52,x53,x54,x56,x61,x62,x63,x64,x65,m123,m124,m125,m126,m134,m135,m136,m145,m146,m156,m234,m235,m236,m245,m246,m256,m345,m346,m356,m456,
MonomialOrder=>Eliminate 30];

M = matrix {
{x12, x13, x14, x15, x16, m123, m124, m125, m126, m134, m135, m136, m145, m146, m156},
{x21, m123, m124, m125, m126, x23, x24, x25, x26, m234, m235, m236, m245, m246, m256},
{m123, x31, m134, m135, m136, x32, m234, m235, m236, x34, x35, x36, m345, m346, m356},
{m124, m134, x41, m145, m146, m234, x42, m245, m246, x43, m345, m346, x45, x46, m456},
{m125, m135, m145, x51, m156, m235, m245, x52, m256, m345, x53, m356, x54, m456, x56},
{m126, m136, m146, m156, x61, m236, m246, m256, x62, m346, m356, x63, m456, x64, x65}};

I = minors(3,M);
time J = ideal selectInSubring(1, gens gb(I,DegreeLimit=>3));
betti mingens J
time J = ideal selectInSubring(1, gens gb(I,DegreeLimit=>5));
betti mingens J
time codim J, degree J



