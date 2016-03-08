import sys
import os
import subprocess

if len(sys.argv) != 2:
    print "Usage: python run_deletion_experiments.py <shell_script>"
    exit(0)


shell_script = sys.argv[1]
#NT = 16
NT = 4
#threadArray = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 35, 40]
threadArray = [0, 1, 2, 3, 4]

os.system("sed -i 's/{[0-9]*}/{0}/g' ../sbcMain.c")
os.system("rm tmp*")
for i in range(1, NT + 1):
    original_thread_string = "{" + str(threadArray[i - 1]) + "}"
    new_thread_string = "{" + str(threadArray[i]) + "}"
    print new_thread_string 
    sed_command = "sed -i 's/" + original_thread_string + "/" + new_thread_string + "/g' ../sbcMain.c"
    os.system(sed_command)
    experiment_command = "./" + shell_script + " > tmp2.csv"
    os.system(experiment_command)
    serialize_command = "grep -vE \"(make|gcc|NV|initial|lm)\" tmp2.csv >> tmp.csv"
    os.system(serialize_command)
