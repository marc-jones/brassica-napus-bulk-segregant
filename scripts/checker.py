import subprocess

md5sums = subprocess.check_output(
    ['find', '-name', 'MD5.txt']).strip().split('\n')

for md5sum in md5sums:
    md5sum_dict = {line.strip().split()[1]: line.strip().split()[0] for line in
        open(md5sum)}
    for key in md5sum_dict.keys():
        actual_md5 = subprocess.check_output(['md5sum',
            md5sum.replace('MD5.txt', key)]).split()[0]
        if actual_md5 == md5sum_dict[key]:
            print('\t'.join([key, 'Success!']))
        else:
            print('\t'.join([key, 'Failure!']))
