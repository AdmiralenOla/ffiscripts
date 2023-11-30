#!/usr/bin/env python

'''This script evaluates taxid membership in contaminants and human-related lists'''

import sys
import pandas as pd


human_related = [24, 158, 164, 192, 196, 199, 200, 203, 204, 205, 206, 238, 239, 247, 253, 263, 287, 292, 293, 303, 306, 316, 329, 358, 391, 446, 447, 449, 450, 453, 459, 463, 470, 472, 477, 478, 479, 480, 483, 485, 486, 487, 488, 489, 490, 495, 498, 502, 504, 511, 518, 519, 520, 529, 539, 540, 545, 550, 562, 565, 569, 571, 573, 577, 582, 585, 587, 614, 615, 616, 618, 630, 644, 646, 648, 650, 652, 663, 666, 669, 670, 672, 674, 676, 678, 714, 719, 723, 726, 727, 729, 730, 732, 735, 739, 740, 747, 749, 752, 753, 754, 760, 813, 817, 824, 837, 850, 851, 860, 861, 885, 1014, 1015, 1017, 1018, 1019, 1034, 1245, 1246, 1260, 1261, 1270, 1273, 1276, 1280, 1282, 1283, 1286, 1288, 1290, 1292, 1296, 1302, 1303, 1304, 1305, 1306, 1308, 1309, 1310, 1311, 1313, 1314, 1318, 1328, 1343, 1348, 1351, 1352, 1353, 1379, 1383, 1396, 1397, 1401, 1402, 1409, 1421, 1465, 1472, 1490, 1498, 1502, 1504, 1506, 1509, 1513, 1529, 1536, 1547, 1559, 1590, 1591, 1596, 1597, 1613, 1622, 1624, 1632, 1639, 1648, 1652, 1655, 1656, 1659, 1660, 1661, 1667, 1685, 1689, 1698, 1701, 1710, 1713, 1717, 1720, 1725, 1747, 1750, 1753, 1766, 1778, 1780, 1781, 1785, 1795, 1824, 1831, 1898, 1931, 2014, 2035, 2047, 2054, 2055, 2104, 2702, 2762, 3074, 5068, 5480, 10279, 10566, 13076, 13689, 28035, 28037, 28038, 28083, 28085, 28087, 28090, 28091, 28095, 28110, 28112, 28119, 28125, 28126, 28129, 28131, 28132, 28133, 28134, 28136, 28137, 28141, 28172, 28188, 28189, 28214, 28264, 28449, 28450, 28454, 29313, 29317, 29321, 29341, 29354, 29363, 29379, 29380, 29382, 29388, 29391, 29394, 29432, 29446, 29466, 29506, 29523, 31312, 31973, 32002, 32013, 33007, 33010, 33028, 33033, 33034, 33053, 33889, 33918, 33935, 33968, 34059, 34062, 34105, 34391, 35519, 35703, 35783, 36740, 36809, 36894, 37326, 37329, 37330, 37637, 37923, 38284, 38286, 38287, 38289, 38290, 38301, 38303, 38304, 38313, 38875, 39777, 39778, 39791, 39950, 40091, 40214, 40216, 40324, 40542, 40543, 41200, 41202, 41276, 41880, 41976, 41986, 42817, 42895, 43263, 43675, 43765, 43767, 43768, 43769, 43770, 43990, 43997, 44058, 44737, 44750, 45242, 45243, 45254, 45634, 46124, 46125, 46771, 47312, 47478, 47671, 47715, 47770, 47846, 47884, 47917, 47920, 48296, 51160, 51671, 52132, 52227, 52768, 52769, 52773, 53363, 53364, 53437, 53462, 53972, 55211, 55565, 56774, 56811, 57171, 58172, 59505, 59561, 59823, 60133, 60552, 61645, 65058, 66228, 66851, 68766, 68891, 68892, 69392, 69823, 69968, 70348, 71254, 71421, 71451, 71999, 72556, 72557, 74319, 74706, 76122, 76123, 76124, 76759, 76773, 76775, 76777, 76832, 76856, 76857, 76859, 77917, 78258, 78259, 78355, 79263, 80866, 80878, 81858, 81950, 82203, 82347, 82380, 82541, 82633, 83558, 84109, 84112, 84163, 84521, 84698, 85698, 86185, 87883, 89093, 90239, 90241, 90245, 93061, 93062, 93218, 93219, 93220, 94138, 95486, 105219, 106648, 106654, 108980, 109790, 110505, 110845, 111015, 113107, 113287, 114527, 114528, 114702, 117179, 119676, 120957, 122586, 122587, 123899, 128780, 131110, 131111, 132249, 132933, 133926, 134034, 134375, 134533, 134534, 134537, 135079, 135080, 135082, 135083, 135487, 137732, 142586, 143361, 143393, 146827, 147645, 150055, 152331, 154288, 156978, 156979, 157687, 157688, 157691, 158822, 158823, 158877, 158878, 158879, 160490, 160491, 161879, 161890, 161902, 169292, 170187, 171101, 176090, 176279, 176280, 177972, 178214, 180332, 180588, 181487, 182337, 184870, 185639, 186103, 187101, 187491, 189423, 189722, 189723, 190304, 192066, 193461, 193567, 194702, 195105, 196620, 197575, 198466, 199591, 200476, 202789, 203275, 206043, 206506, 207340, 209882, 210007, 216816, 217203, 217204, 218538, 220685, 220687, 221027, 223392, 225324, 228599, 228603, 228604, 230143, 237576, 240125, 242619, 243275, 243701, 244292, 246198, 249188, 257758, 262727, 262728, 267212, 267747, 272548, 272556, 272622, 272623, 272831, 273036, 273136, 280147, 281310, 281920, 282305, 282402, 282458, 282459, 283734, 285091, 286636, 286802, 291112, 293653, 305719, 310300, 319701, 319706, 319709, 319939, 321967, 322095, 324831, 326522, 326523, 327574, 327575, 334390, 341694, 341722, 342002, 347253, 352165, 354243, 354351, 359786, 359787, 360104, 360105, 361500, 362948, 363952, 370551, 370552, 370553, 370554, 370622, 371601, 373153, 374833, 374927, 374928, 374930, 374931, 374932, 374933, 375063, 375177, 375432, 379413, 386414, 388919, 393480, 399795, 400946, 406556, 406557, 406558, 406559, 406560, 406561, 406562, 406563, 407975, 411465, 411466, 411570, 411577, 416870, 418127, 419015, 419208, 419475, 420404, 426430, 431269, 431947, 435830, 435832, 435838, 439703, 447455, 447456, 450394, 451515, 451516, 452948, 453361, 453362, 453363, 453364, 453365, 453366, 455227, 456482, 457403, 457405, 457921, 461393, 467210, 467705, 469378, 469599, 469601, 469602, 469604, 469607, 469621, 470565, 471872, 471876, 472693, 479117, 479436, 480035, 486408, 487213, 487214, 487215, 488221, 488222, 488223, 488730, 489653, 491076, 493803, 497962, 497963, 497980, 501496, 502790, 511691, 512566, 512767, 512768, 512769, 516950, 520603, 521004, 521005, 521095, 521097, 521392, 521393, 521520, 523794, 525283, 525325, 525326, 525337, 525361, 525374, 525375, 525376, 525378, 525381, 537973, 544580, 544581, 546262, 546263, 546264, 546265, 546266, 546268, 546269, 546270, 546271, 546273, 546274, 546275, 546342, 547045, 548470, 548473, 548474, 548475, 553171, 553174, 553175, 553177, 553178, 553184, 553198, 553199, 553201, 553207, 553218, 553219, 553220, 553565, 553567, 553568, 553571, 553573, 553574, 553577, 553580, 553581, 553583, 553588, 553590, 553592, 553594, 553596, 553601, 554406, 556263, 556499, 561276, 562973, 562981, 562982, 562983, 563032, 563033, 566549, 568703, 568704, 574093, 575590, 575593, 575611, 575612, 575614, 575615, 585161, 585501, 585503, 589436, 591365, 592010, 592026, 592028, 592031, 595501, 596085, 596315, 596317, 596319, 596320, 596322, 596323, 596324, 596329, 596330, 604162, 607712, 608534, 619693, 620833, 620903, 626084, 626369, 626522, 626523, 629741, 630527, 630588, 633147, 633701, 634176, 634994, 638300, 638301, 638302, 638849, 641143, 641147, 641149, 645512, 649742, 649743, 649760, 649761, 649764, 651822, 652722, 653386, 655813, 656912, 656913, 662598, 665914, 668336, 671211, 671214, 671218, 671224, 678932, 679188, 679192, 679193, 679194, 679195, 679196, 679198, 679199, 679200, 679201, 680646, 681288, 684066, 684738, 686659, 686660, 693991, 694569, 696216, 699187, 702437, 702438, 702439, 702745, 703339, 706433, 706434, 706436, 706437, 706439, 712116, 712122, 712150, 712310, 712357, 712361, 712362, 712363, 712365, 712368, 712411, 712435, 712466, 712471, 712528, 712538, 712623, 712624, 712633, 712710, 712711, 712938, 712961, 712982, 712991, 713030, 713051, 713059, 714315, 742814, 742820, 746361, 746830, 746831, 746832, 748671, 749551, 754505, 754506, 754507, 755171, 755172, 756689, 759851, 760570, 760746, 760787, 760791, 760809, 760810, 760834, 760861, 762948, 762963, 762965, 764544, 767029, 767031, 767100, 767453, 768724, 768726, 768727, 768728, 796937, 796942, 796943, 797473, 798300, 857099, 857100, 857101, 857102, 857103, 857104, 857105, 857106, 857107, 857108, 857109, 857110, 857111, 857112, 857113, 857114, 857115, 857116, 857117, 857118, 857119, 857120, 857121, 857122, 857123, 857124, 857125, 857126, 857127, 857128, 857129, 857130, 857131, 857132, 857133, 857134, 857135, 857136, 857137, 857138, 857139, 857140, 857141, 857142, 857143, 857144, 857145, 857146, 857147, 857148, 857149, 857150, 857151, 857152, 857153, 857154, 857155, 857290, 857291, 857292, 857571, 857572, 857573, 857574, 857575, 857576, 857577, 857578, 857579, 857581, 861450, 861452, 861454, 861455, 862513, 862515, 862964, 862966, 862967, 862968, 862969, 862970, 862971, 864563, 864565, 864567, 864568, 864570, 866630, 866776, 866778, 868129, 869214, 869215, 869216, 869269, 869309, 871237, 871541, 873513, 873517, 873533, 879309, 879310, 880592, 883092, 883094, 883103, 883109, 883158, 883167, 885272, 886289, 887325, 887898, 887901, 887929, 888019, 888048, 888049, 888050, 888051, 888052, 888054, 888055, 888056, 888057, 888059, 888060, 888061, 888062, 888721, 888727, 888728, 888741, 888742, 888743, 888746, 888808, 888809, 888810, 888811, 888812, 888813, 888814, 888815, 888816, 888825, 888832, 888833, 889201, 889204, 889206, 904294, 904296, 904306, 904317, 904338, 905067, 907486, 907487, 907488, 907489, 907490, 907491, 907492, 907493, 908937, 909420, 909952, 912594, 927666, 929102, 929793, 935589, 935598, 935599, 935897, 936561, 936563, 936589, 936596, 942513, 944557, 944560, 944564, 944565, 945844, 947033, 947828, 979627, 985002, 985008, 997347, 997352, 997353, 997356, 997830, 999414, 999415, 999422, 999423, 999424, 999425, 999426, 999427, 999428, 999429, 999430, 999431, 999432, 999433, 999434, 999435, 999436, 999437, 999438, 999439, 999440, 1000570, 1000588, 1000590, 1002365, 1005704, 1005705, 1009852, 1028802, 1028803, 1028804, 1028805, 1028806, 1029822, 1030843, 1031709, 1032505, 1035184, 1035185, 1035188, 1035189, 1035190, 1035193, 1035194, 1035195, 1035196, 1035197, 1041521, 1042402, 1046624, 1046629, 1047171, 1048332, 1051006, 1051972, 1051985, 1069623, 1069625, 1069626, 1069628, 1073353, 1073362, 1073366, 1073367, 1073372, 1074066, 1074092, 1074093, 1074095, 1074100, 1074101, 1074102, 1074104, 1074105, 1074106, 1074107, 1074108, 1074109, 1074111, 1074112, 1074113, 1074114, 1074115, 1074116, 1074118, 1074119, 1074120, 1074121, 1074122, 1074123, 1074124, 1074125, 1074126, 1074127, 1074128, 1074129, 1074130, 1074134, 1074135, 1074136, 1074137, 1074138, 1074140, 1074143, 1074144, 1074146, 1074148, 1074149, 1074151, 1074153, 1074154, 1074155, 1074156, 1074157, 1074159, 1074160, 1074161, 1074162, 1074163, 1074164, 1074165, 1074166, 1074167, 1074168, 1074169, 1074170, 1074171, 1074173, 1074175, 1074176, 1074177, 1074178, 1074179, 1074180, 1074181, 1074182, 1074183, 1074184, 1074185, 1074186, 1074190, 1074494, 1078480, 1078483, 1081904, 1088720, 1089447, 1091045, 1095729, 1095730, 1095731, 1095733, 1095738, 1095739, 1095740, 1095741, 1095742, 1095743, 1095744, 1095747, 1095748, 1095750, 1095752, 1104322, 1108963, 1110546, 1111678, 1114965, 1114966, 1114967, 1114969, 1115803, 1115809, 1120941, 1120942, 1120943, 1120944, 1120957, 1120979, 1121268, 1121367, 1122171, 1122172, 1122174, 1122949, 1122980, 1122982, 1122984, 1122985, 1122986, 1122987, 1122989, 1122991, 1122993, 1122994, 1123249, 1123263, 1123310, 1123317, 1125699, 1125700, 1125701, 1125702, 1125712, 1125717, 1125718, 1125719, 1125720, 1125721, 1125722, 1125723, 1125724, 1125725, 1127690, 1127691, 1127692, 1127693, 1127694, 1127695, 1127696, 1127699, 1128111, 1130804, 1138874, 1141657, 1155071, 1157946, 1159208, 1161421, 1161422, 1161424, 1161902, 1167007, 1167008, 1167009, 1167010, 1167628, 1177574, 1177728, 1185324, 1190621, 1194526, 1195243, 1198676, 1200793, 1203258, 1203259, 1203550, 1203557, 1203559, 1203561, 1203562, 1203566, 1203602, 1203603, 1203619, 1203622, 1203624, 1203625, 1203627, 1203632, 1211023, 1216362, 1225186, 1225187, 1225188, 1225189, 1225190, 1225191, 1225192, 1225193, 1225194, 1225195, 1225196, 1225197, 1225198, 1225199, 1225200, 1225201, 1225202, 1225203, 1225204, 1225205, 1226633, 1227261, 1227262, 1227264, 1227266, 1227268, 1227269, 1227270, 1227271, 1227272, 1227276, 1234601, 1234680, 1234877, 1235815, 1236497, 1236504, 1236508, 1236516, 1236517, 1236518, 1236608, 1239307, 1241978, 1242967, 1242969, 1243032, 1244083, 1248420, 1256219, 1256223, 1256230, 1257037, 1257038, 1257039, 1257040, 1257041, 1257042, 1266996, 1266997, 1273133, 1283280, 1287474, 1287476, 1287736, 1292047, 1292048, 1293577, 1297564, 1297565, 1297566, 1297567, 1302863, 1307427, 1307428, 1307442, 1307443, 1307444, 1309795, 1311575, 1316254, 1316593, 1316596, 1316933, 1318634, 1321772, 1321774, 1321775, 1321779, 1321781, 1321782, 1321784, 1321786, 1321814, 1321815, 1321816, 1321817, 1321818, 1321820, 1321821, 1321822, 1321823, 1334627, 1346615, 1347368, 1347369, 1347790, 1353243, 1365628, 1366052, 1380685, 1389713, 1389922, 1395125, 1401068, 1401072, 1401073, 1401077, 1401079, 1403335, 1403336, 1403338, 1403829, 1403949, 1404260, 1407647, 1410950, 1411021, 1411022, 1411148, 1411915, 1415626, 1423782, 1423799, 1423814, 1430326, 1434258, 1434259, 1434260, 1434261, 1434262, 1434263, 1434264, 1434265, 1437447, 1440768, 1440769, 1440770, 1440771, 1448849, 1448850, 1501329, 1501332, 1522312, 1544413, 1544416, 1547448, 1578165, 1583331, 1661745, 1671022, 1671023, 1673725, 1697053, 1705617, 1715123, 1715217, 1739279, 1739317, 1739435, 1739543, 1756149, 1776741, 1785995, 1817405, 1834153, 1852361, 1859694, 1869190, 1871047, 1871052, 1872515, 1874826, 1884263, 1889813, 1890675, 1891233, 1891644, 1895474, 1911679, 1924944, 1944660, 1960874, 1965292, 1977869, 1986155, 2079439, 2081702, 2081962, 2093824, 2094119, 2126346, 2382124, 2510778, 2572087, 2572088, 2572089, 2748177, 2748316, 2748317]


contaminants = [10, 16, 20, 22, 68, 75, 81, 88, 157, 171, 194, 222, 237, 239, 245, 247, 258, 265, 283, 285, 286, 293, 294, 296, 303, 304, 316, 329, 338, 357, 358, 374, 375, 376, 379, 382, 407, 410, 434, 469, 470, 471, 482, 483, 484, 497, 505, 507, 528, 532, 547, 561, 562, 564, 570, 573, 613, 615, 623, 713, 724, 726, 729, 816, 818, 820, 821, 823, 836, 838, 840, 846, 851, 853, 906, 958, 963, 986, 991, 1016, 1017, 1033, 1076, 1243, 1257, 1269, 1270, 1279, 1282, 1283, 1290, 1298, 1301, 1302, 1303, 1304, 1305, 1308, 1309, 1314, 1338, 1343, 1350, 1351, 1352, 1357, 1358, 1375, 1379, 1380, 1382, 1386, 1402, 1409, 1485, 1492, 1522, 1547, 1567, 1578, 1613, 1624, 1654, 1655, 1656, 1660, 1663, 1678, 1696, 1698, 1716, 1743, 1747, 1785, 1827, 1828, 1833, 1835, 1839, 1847, 1860, 2034, 2040, 2060, 2702, 2717, 2718, 2736, 2741, 2742, 2745, 2755, 5308, 5665, 9606, 10264, 10288, 10298, 10345, 10359, 10376, 10583, 10665, 10732, 10756, 10760, 10761, 10798, 10847, 10863, 10864, 11033, 11103, 11620, 11676, 11801, 11809, 11861, 11864, 11866, 11867, 11870, 11874, 11876, 11878, 11884, 11885, 11886, 11888, 11908, 11946, 11948, 11950, 11955, 11958, 11960, 11987, 12056, 12071, 12104, 12138, 12235, 12916, 12960, 13275, 13687, 28035, 28037, 28050, 28090, 28100, 28101, 28116, 28117, 28124, 28127, 28132, 28135, 28214, 28344, 28453, 28901, 29330, 29363, 29389, 29404, 29448, 29465, 29466, 29536, 29570, 29580, 29581, 31535, 31552, 31647, 31668, 31669, 31670, 31988, 31998, 32008, 32067, 32199, 32207, 32257, 33010, 33011, 33029, 33033, 33038, 33042, 33050, 33057, 33747, 33870, 33882, 34004, 34062, 34072, 34073, 35306, 35812, 35841, 36381, 36773, 37914, 38303, 38304, 38305, 38313, 39486, 39491, 39492, 39643, 39778, 39791, 39948, 40214, 40215, 40323, 40324, 40520, 40542, 40979, 41275, 41276, 43668, 43768, 43992, 43994, 44249, 45634, 45669, 46123, 46124, 46125, 46205, 46353, 46466, 46506, 46913, 47494, 47496, 47883, 47920, 48736, 50340, 50709, 51101, 51680, 52133, 52972, 53370, 53412, 53457, 54007, 54914, 55080, 55197, 55508, 55510, 56946, 57493, 57495, 59732, 59753, 59803, 60550, 61654, 64001, 64974, 66831, 68287, 68347, 68569, 68887, 68892, 68909, 69359, 69392, 69823, 70774, 70863, 71999, 72000, 72201, 74030, 74316, 74426, 75654, 75659, 76731, 76761, 76773, 76775, 76890, 77097, 77583, 80865, 80866, 80878, 80882, 82380, 83618, 84108, 84112, 84135, 84292, 84567, 84756, 85413, 85506, 85698, 85708, 86182, 86669, 87883, 89966, 90963, 92442, 92793, 93064, 93681, 94008, 94625, 94626, 96345, 97050, 98513, 99158, 100175, 102148, 106588, 106589, 106592, 106648, 106649, 108981, 114248, 114702, 115553, 115555, 115979, 115987, 117207, 117506, 117563, 119852, 120831, 125216, 129337, 129817, 129951, 131079, 134533, 134534, 135517, 135858, 136084, 136273, 136996, 137732, 138336, 139872, 144191, 145579, 146500, 146937, 147207, 147645, 149698, 150247, 151416, 154046, 154288, 155892, 157076, 159150, 160404, 161492, 161879, 162426, 165179, 165695, 165696, 165697, 165779, 169292, 169683, 169864, 172042, 172044, 172088, 172371, 172827, 174708, 176652, 180282, 181522, 182269, 185636, 185638, 185950, 186617, 188350, 190721, 192843, 194802, 198112, 198252, 201096, 201847, 202952, 202954, 202956, 204456, 204516, 204525, 205844, 205879, 208479, 211589, 212035, 212791, 213484, 215579, 215796, 216465, 216816, 216851, 218538, 220638, 222991, 225324, 225991, 227470, 228654, 230120, 230123, 231454, 232523, 233894, 236752, 238854, 239935, 242521, 242527, 242861, 244127, 246602, 246787, 249058, 251749, 255204, 256325, 257758, 260149, 260373, 261299, 263375, 266749, 269069, 269447, 269448, 270673, 271647, 272239, 274591, 278008, 279280, 281915, 285986, 287412, 293256, 294382, 294631, 295418, 299566, 301302, 308865, 309120, 310297, 310298, 320843, 320850, 321895, 322019, 323621, 328552, 331278, 332102, 333297, 334852, 335058, 335924, 336724, 337315, 340016, 343873, 345198, 346179, 346932, 347326, 347327, 352475, 354090, 354354, 356778, 358705, 361607, 362076, 362413, 362418, 363020, 363745, 369581, 369926, 370959, 370974, 373126, 374425, 376469, 376758, 378210, 381630, 386793, 387661, 399781, 400567, 401469, 401472, 414970, 418240, 423604, 428988, 431058, 431059, 432308, 432330, 435913, 437755, 440250, 444860, 444861, 444862, 445688, 445700, 446529, 455370, 462590, 467210, 469322, 469660, 475299, 490913, 493803, 501783, 503361, 504481, 504501, 504553, 518981, 529883, 529884, 536444, 536454, 536473, 537874, 541865, 551895, 554168, 555387, 560405, 568987, 572511, 573173, 573174, 573176, 577310, 582345, 585044, 590739, 642253, 642255, 644524, 664785, 665874, 669464, 682522, 682650, 691965, 693272, 697906, 742919, 743583, 743813, 745100, 745102, 745107, 747294, 749413, 753084, 753085, 754042, 754044, 754048, 754051, 754052, 754058, 754060, 754067, 754072, 754075, 756275, 756277, 756279, 756280, 756282, 756892, 757342, 759804, 760732, 763552, 764562, 765765, 863372, 866889, 877240, 880162, 889876, 889949, 889956, 908819, 925983, 925984, 926067, 926697, 929832, 929835, 938080, 938081, 938082, 944322, 946234, 948071, 948870, 979525, 979534, 981335, 987053, 989370, 990721, 993502, 994601, 998086, 1004300, 1004302, 1010698, 1026955, 1034111, 1034128, 1034806, 1036779, 1042123, 1043493, 1045778, 1048515, 1048516, 1048517, 1048520, 1048521, 1054211, 1054461, 1055192, 1056830, 1072204, 1074214, 1075768, 1076759, 1084719, 1089119, 1094892, 1112209, 1114179, 1127514, 1127515, 1128131, 1128140, 1128143, 1128151, 1128422, 1129146, 1129191, 1131316, 1131812, 1132026, 1133022, 1133292, 1134405, 1137745, 1141134, 1141135, 1141136, 1147094, 1147148, 1148801, 1150298, 1150989, 1150991, 1154691, 1159870, 1160721, 1161935, 1162295, 1168281, 1168478, 1168479, 1168549, 1169627, 1170653, 1173749, 1173759, 1173761, 1173762, 1176422, 1176423, 1187128, 1188795, 1195080, 1197951, 1198452, 1204514, 1204517, 1204529, 1204539, 1208587, 1211417, 1217692, 1222338, 1224510, 1225745, 1229753, 1229760, 1229782, 1229784, 1229786, 1229787, 1229788, 1229789, 1229790, 1229791, 1229792, 1229793, 1229794, 1230476, 1231048, 1235314, 1235647, 1235648, 1235649, 1235650, 1235653, 1235654, 1235655, 1235656, 1235657, 1235689, 1237364, 1243183, 1262072, 1262517, 1264700, 1269028, 1276755, 1278246, 1278247, 1278248, 1278249, 1278250, 1278251, 1278252, 1278255, 1278261, 1278263, 1278265, 1278278, 1279082, 1283071, 1283076, 1283077, 1285382, 1296654, 1299317, 1301280, 1316739, 1320556, 1325731, 1327037, 1327964, 1327970, 1327971, 1327972, 1327974, 1327975, 1327976, 1327977, 1327979, 1327980, 1327981, 1327982, 1327983, 1327985, 1327990, 1327992, 1327993, 1327995, 1328029, 1328030, 1332312, 1335230, 1340709, 1340810, 1340812, 1341019, 1348912, 1352534, 1353941, 1357423, 1357714, 1379694, 1379788, 1385658, 1385659, 1391223, 1394983, 1395610, 1395611, 1395612, 1395613, 1395614, 1395615, 1395616, 1395617, 1395618, 1395619, 1395620, 1399915, 1406795, 1407671, 1414739, 1414766, 1416009, 1416628, 1416631, 1417599, 1420594, 1424633, 1429767, 1434319, 1435036, 1435438, 1449437, 1450746, 1450749, 1454023, 1454024, 1454025, 1455074, 1458842, 1458843, 1458858, 1458859, 1458860, 1458861, 1461100, 1461743, 1462581, 1465639, 1472912, 1474867, 1476390, 1477406, 1478972, 1481186, 1481785, 1486472, 1486657, 1486662, 1492737, 1492738, 1497615, 1497851, 1498188, 1500757, 1504823, 1505225, 1506553, 1509403, 1511763, 1519385, 1519387, 1519389, 1519390, 1519395, 1519396, 1519397, 1519405, 1519439, 1521385, 1521387, 1521389, 1524880, 1524881, 1524882, 1526550, 1527506, 1527515, 1527519, 1527524, 1529058, 1536592, 1537091, 1538644, 1540094, 1540097, 1540098, 1540099, 1541891, 1548900, 1548901, 1548918, 1557033, 1560342, 1566990, 1567016, 1567475, 1574422, 1589298, 1592081, 1592083, 1592085, 1592086, 1592088, 1592093, 1592095, 1592096, 1592112, 1592113, 1592127, 1592207, 1592212, 1608440, 1608451, 1609634, 1618248, 1618254, 1619232, 1646498, 1654716, 1665556, 1678129, 1678143, 1706231, 1714344, 1736272, 1736280, 1736296, 1736316, 1736528, 1736532, 1740090, 1835254, 1843761, 1891767, 1895944, 1977402, 2282523, 2559073]

with open(sys.argv[1],'r') as my_infile:
     my_infile_pd = pd.read_csv(my_infile, sep="\t")
     print(my_infile_pd)
     my_infile_pd["Is_contaminant"] = my_infile_pd["taxonomy_id"].apply(lambda x: 1 if x in contaminants else 0)
     my_infile_pd["Is_human_related"] = my_infile_pd["taxonomy_id"].apply(lambda x: 1 if x in human_related else 0)
     outfile = sys.argv[1].rstrip(".tsv") + "_.tsv"
     with open(outfile,'w') as my_outfile:
        my_infile_pd.to_csv(my_outfile,sep="\t",index=False)