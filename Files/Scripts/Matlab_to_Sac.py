import scipy
import os

file_path = "/Users/stephenyoung/EQ_Stuff/EQ_Files/Data_sets/06_25_2019"

file_list = os.listdir(file_path)
print(scipy.io.whosmat('/Users/stephenyoung/EQ_Stuff/EQ_Files/Data_sets/06_25_2019/06252019_154340_UTC.mat'))
june_six = scipy.io.loadmat('/Users/stephenyoung/EQ_Stuff/EQ_Files/Data_sets/06_25_2019/06252019_154340_UTC.mat')
print(june_six)

# ['data-ready-for-viewing']



# for file_name in file_list:
#     file_path = os.path.join(file_name)
#     if os.path.isfile(file_path):
#         scipy.io.loadmat(file_path)





