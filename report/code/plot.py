from matplotlib import pyplot as plt
from math import sqrt

x = range(0, 256)

yr = [118730, 124423, 144221, 1225318, 340183, 1213479, -590360, -3497253, -120686, -94267, -6444, -804154, -236760, -2438917, 416198, -542322, -38733, -27901, -18985, -96706, -7221, 242768, -53983, -261752, -17173, -102787, -317369, 1080021, 287987, 143900, 111449, 1620936, -26336, -22865, -19619, 39576, -13883, -105748, 19835, 2466, -5811, 11152, 47576, -158148, -50364, -127484, 14888, -240149, -8119, -12151, -16319, -47807, -25305, -510213, 49399, -84278, -55436, -77147, -124744, 252678, 13977, 461034, -93969, -25990, -23548, -22472, -21392, 89740, -19242, -199603, 29344, 96728, -14992, 3908, 31957, -102502, -20034, 424632, -108972, -60610, -6756, -5391, -3391, 106128, 3036, 257706, 1602, 160721, 38263, 75029, 176099, -568087, -131473, -337765, -7460, -719724, 7781, 7818, 7788, -29490, 7553, 118073, -504, -17067, 6621, 12384, 35804, -85525, -35230, -89583, 3584, -220955, 4606, 3523, 2447, -41437, 515, 43175, 3963, -37138, 871, -9679, -22292, -19189, -9764, 15454, -24004, -43032, -22990, -22991, -22993, 73816, -22997, -242913, 10776, 26701, -22992, -21841, -33964, 16720, 4648, 296133, -68534, -20306, -23013, -20819, -18207, 79462, -11445, -37688, 55095, 113114, 16595, 37577, 105415, -456267, -110845, -60092, -54397, -394808, -23086, -22263, -21225, 113012, -18503, -166754, 39511, 158598, -10111, 26146, 103046, -261954, -69324, 625178, -187884, -366081, 23619, 27425, 31243, -15789, 39825, 966225, -67125, 36786, 81266, 81859, 111602, -371502, -40939, -530502, 55639, -76828, -23548, -24632, -25704, 55966, -27900, -290519, -3644, -39632, -32214, -42608, -79897, 72264, 14082, 246048, -57556, -17454, -40742, -37453, -33497, 108452, -23030, -141958, 107538, 172605, 24119, 58405, 181179, -889931, -225593, -66515, -139056, -922256, -57203, -61578, -65884, 160902, -74411, -845183, 11562, -14997, -91231, -101898, -174870, 71643, 45134, 1325001, -337228, 341773, -129350, -146273, -163859, 443513, -204103, -2952383, 195359, 216046, -379921, -473337, -849034, 4184461, 507570, 5189718, -392378, 1974974]
yi = [0, 72632, 164733, -3130100, 716155, 2189719, -323828, 1504482, -386260, -372249, -476822, 673942, -53957, 2728810, -567228, -1157610, -125797, -110685, -94576, 376133, -58317, -164107, 195095, 449224, 61220, 129271, 333785, -1264421, -302876, -307301, -95477, -962980, -55652, -52205, -48486, 241037, -40260, -387174, 80407, 304090, -20649, 45868, 178900, -448001, -115819, 1042572, -304127, -601930, 34581, 39202, 43534, -14332, 52463, 1237108, -87672, 46535, 94075, 91085, 115400, -345002, -31029, -516246, 59137, -16119, -22610, -22033, -21420, 111543, -20060, -199802, 30548, 132769, -16863, 1199, 24933, -63469, -11521, 635373, -155521, -60950, -8246, -7179, -5308, 143284, 1844, 241518, 19937, 255414, 51065, 118720, 307534, -1047112, -262056, -492812, -50428, -1600255, 21535, 22821, 23986, -26191, 25989, 419034, -13817, 20421, 28874, 25739, 41898, 55347, -23664, -218459, 58245, -352670, 32572, 34187, 35221, -59817, 35933, 414897, -37376, -25699, 35754, 38984, 51768, -177454, -3406, -175621, 18986, 71395, 0, 280, 553, 30832, 1107, 14505, 17706, 81706, 2218, 15333, 40536, -62054, -11731, 267860, -67648, -29968, 4481, 3767, 3440, 53575, 4315, 126581, -1953, 101936, 20062, 48389, 125451, -414571, -105674, -176651, -24479, -652138, 9334, 8365, 7428, -3735, 5686, 114256, -14671, 3404, 3317, -7534, -26806, 119579, 18701, 128624, -20587, -27140, 9489, 12466, 15554, 36982, 22361, 406526, -46426, 63301, 46267, 65785, 105288, -205258, -4545, -454830, 96851, 58919, 22610, 23153, 23658, -44475, 24558, 289246, 1110, 32015, 25941, 34001, 66517, -63839, -21391, -163175, 33221, -82504, 27178, 23989, 20592, -59592, 13094, 113496, -66727, -84070, -10071, -27164, -88994, 464852, 112962, -22648, 80856, 421561, 24783, 25715, 26500, -54871, 27665, 303656, -5975, 9233, 28526, 27659, 39916, 11895, -2850, -315989, 80837, -35860, 25742, 27349, 28591, -66225, 30395, 389005, -24118, -20953, 36524, 40626, 62776, -251778, -24072, -184131, 9042, -23639]

y = [sqrt(r**2 + i**2) for r, i in zip(yr, yi)]

ax = plt.axes(projection = '3d')
ax.set_proj_type('ortho')

ax.plot3D(x, yr, yi)
plt.show()
