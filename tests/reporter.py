# #! python3
from bs4 import BeautifulSoup
import os,sys #.path import path
from colorama import Fore as colors

libs = {
    'bloch':['bloch_base','continuous','felix'],
    'EDutils':['pets']
}
libs = {
    'bloch':['bloch_base','continuous','felix'],
    'EDutils':['pets']
}

def file_summary(file):
    body=open(file,'rb')
    soup = BeautifulSoup(body,'html.parser')

    # summary = "- [%s](%s) : " %(os.path.basename(file),reports_lnk(file))
    # summary+=  ','.join([ soup.find_all('span',s)[0].get_text()
    #     for s in ['passed','skipped','failed']])

    summary = "[%s](%s) | " %(os.path.basename(file),reports_lnk(file))
    summary+=  '|'.join([ soup.find_all('span',s)[0].get_text().split(' ')[0]
        for s in ['passed','skipped','failed']])

    return summary

def write_report():
    report="""
    # Report
    """

    for lib,files in libs.items():
        report+="\n## %s\n" %lib
        report+="test | passed | skipped | failed  \n-- | -- | -- | -- \n"
        report+='\n'.join([file_summary('%s/report/test_%s.py.html' %(lib,file))
            for file in files])
        report+="\n"

    report_file="report.md"
    with open(report_file,"w") as f: f.write(report)
    print(colors.YELLOW+"file://"+os.path.realpath(report_file)+colors.RESET)


def run_tests():


if __name__=="__main__":
    hostname=check_output('hostname -A', shell=1).decode().strip().split()[-1]
    # file='EDutils/report/test_pets.py.html'
    # reports_lnk=lambda file:"file://"+os.path.realpath(file)
    # reports_lnk=lambda file:"http://rcccp4s004.rc-harwell.ac.uk:8010/%s" %file
    reports_lnk=lambda file:"http://%s:8010/%s" %(hostname,file)

    run_tests()
    write_report()
