from collections import defaultdict

interActions1 = [

('miR-223', 'ICAM1', 'inhibit'),
('miR-17-3p', 'ICAM1', 'inhibit'),
('miR-31', 'E-selectin', 'inhibit'),
('miR-126', 'VCAM1', 'inhibit'),
('miR-126-5p', 'VCAM1', 'inhibit'),
('miR-126', 'DLK1', 'inhibit'),
('miR-126-5p', 'DLK1', 'inhibit'),
('miR-143', 'JAMA', 'inhibit'),
('miR-145', 'JAMA', 'inhibit'),
('miR-217', 'SIRT1', 'inhibit'),
('miR-155', 'eNOS', 'inhibit'),
('miR-155', 'ETS1', 'inhibit'),
('miR-21', 'PPARA', 'inhibit'),
('miR-712', 'TIMP3', 'inhibit'),
('miR-10a', 'MAP3K7', 'inhibit'),
('miR-181b', 'Importin-a3', 'inhibit'),
('miR-146a', 'TRAF6', 'inhibit'),
('miR-146b', 'TRAF6', 'inhibit'),
('miR-146a', 'IRAK1', 'inhibit'),
('miR-146b', 'IRAK1', 'inhibit'),
('miR-146a', 'IRAK2', 'inhibit'),
('miR-146b', 'IRAK2', 'inhibit'),

('miR-143', 'KLF2', 'activate', 'rev'),
('miR-145', 'KLF2', 'activate','rev'),
('miR-126', 'KLF2', 'activate', 'rev'),
('miR-126-5p', 'KLF2', 'activate', 'rev'),

('miR-126', 'CXCL12', 'activate', 'rev'),
('miR-126-5p', 'CXCL12', 'activate', 'rev'),
('miR-92a', 'KLF2', 'inhibit', 'rev'),
('miR-92a', 'KLF4', 'inhibit', 'rev'),
('miR-92a', 'SOCS5', 'inhibit', 'rev'),
('miR-663', 'KLF4', 'inhibit', 'rev'),

('miR-126', 'ETS1', 'activate', 'rev'),
('miR-126-5p', 'ETS1', 'activate', 'rev'),

]

interActions2 = [

    ('miR-29', 'IFNG', 'inhibit'),
    ('miR-155', 'IFNG', 'inhibit'),
    ('miR-146a', 'CD40', 'inhibit'),
    ('miR-146a', 'CD40L', 'inhibit'),
    ('miR-181a', 'CD40', 'inhibit'),
    ('miR-181a', 'CD40L', 'inhibit'),
    ('miR-21', 'MMP2', 'inhibit'),
('miR-21', 'MMP2', 'inhibit'),
('miR-21', 'MMP3', 'inhibit'),
('miR-21', 'MMP8', 'inhibit'),
('miR-21', 'MMP9', 'inhibit'),
('miR-21', 'MMP12', 'inhibit'),
('miR-21', 'MMP13', 'inhibit'),
('miR-21', 'MMP14', 'inhibit'),

('miR-24', 'MMP2', 'inhibit'),
('miR-24', 'MMP3', 'inhibit'),
('miR-24', 'MMP8', 'inhibit'),
('miR-24', 'MMP9', 'inhibit'),
('miR-24', 'MMP12', 'inhibit'),
('miR-24', 'MMP13', 'inhibit'),
('miR-24', 'MMP14', 'inhibit'),

('miR-29', 'MMP2', 'inhibit'),
('miR-29', 'MMP3', 'inhibit'),
('miR-29', 'MMP8', 'inhibit'),
('miR-29', 'MMP9', 'inhibit'),
('miR-29', 'MMP12', 'inhibit'),
('miR-29', 'MMP13', 'inhibit'),
('miR-29', 'MMP14', 'inhibit'),


]

interactions3 = [
    ('miR-30c', 'MTTP', ''),
('miR-33', 'ABCA1', ''),
('miR-33', 'ABCG1', ''),
('miR-33', 'CPT1A', ''),
('miR-92a', 'KLF2', ''),
('miR-92a', 'KLF4', ''),
('miR-92a', 'SOCS5', ''),
('miR-126-5p', 'DLK1', ''),
('miR-126-3p', 'RGS16', ''),

('miR-145', 'KLF4', ''),
('miR-155', 'BCL6', ''),
('miR-155', 'SOCS1', ''),
('miR-155', 'MAP3K10', ''),

('miR-181b', 'KPNA4', ''),
('miR-302a', 'ABCA1', ''),
('miR-342-5p', 'AKT1', ''),
('miR-712', 'TIMP3', ''),


]


interActions1Dict = defaultdict(list)
for x in interActions1:
    interActions1Dict[x[1]].append(x[0])


interActions2Dict = defaultdict(list)
for x in interActions2:
    interActions2Dict[x[1]].append(x[0])

interActions3Dict = defaultdict(list)
for x in interactions3:
    interActions3Dict[x[1]].append(x[0])

print("andreouInts1=", interActions1Dict)
print("andreouInts2=", interActions2Dict)
print("andreouInts3=", interActions3Dict)



