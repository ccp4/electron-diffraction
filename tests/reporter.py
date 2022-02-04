#!python3
from bs4 import BeautifulSoup
import os,sys #.path import path
from subprocess import check_output,Popen,PIPE
from colorama import Fore as colors

usage='''Usage:
report.py [opt] [libs]

OPTIONS
-------
    opt  : h(help) s(skip running tests) l(libs) f(file serv)
    libs : <lib1>,<lib2>,.. if 'l' in opts
'''


def get_report_summary(file):
    body=open(file,'rb')
    soup = BeautifulSoup(body,'html.parser')
    tests =  '|'.join([ soup.find_all('span',s)[0].get_text().split(' ')[0]
        for s in ['passed','skipped','failed','error']])
    # link_name = os.path.basename(file)#'report'
    summary= "[%s](%s) | %s " %('tests_report',reports_lnk(file),tests)

    return summary

def get_coverage_summary(file):
    body=open(file,'rb')
    soup = BeautifulSoup(body,'html.parser')
    # cov=soup.find_all('tfoot')[0].find('td','right').get_text()
    cov=soup.find_all('span','pc_cov')[0].get_text()
    return " [%s](%s) | %s " %('cov_report',reports_lnk(file),cov)

def run_tests():

    report="# Report\n\n"
    report+="lib | report | passed | skipped | failed | errors | cov_report | coverage \n"
    report+="--  | --     | --     | --      | --     | --     | --         | --       \n"

    for lib in libs:
        report_file = '%s/%s/report/report.html'  %(tests_dir,lib)
        cov_file    = '%s/%s/htmlcov/index.html'  %(tests_dir,lib)

        #run test
        args="--html=%s --self-contained-html --cov-report=html --cov=../../%s" %(report_file,lib)
        cmd ='python3 -m pytest %s' %args
        job ='cd %s/%s; if [ ! -d report ];then mkdir report;fi; %s ' %(tests_dir,lib,cmd)
        print(cmd)
        if run_opt:
            # oe=check_output(job, shell=1).decode().strip()
            p=Popen(job, shell=1);p.wait()
            print(p.communicate())

        #get summary
        lib_report=get_report_summary(report_file)
        lib_cov   =get_coverage_summary(cov_file)

        report+="%s | %s | %s \n" %(lib,lib_report,lib_cov)

    report_file=os.path.join(tests_dir,"report.md")
    with open(report_file,"w") as f: f.write(report)
    print(colors.YELLOW+reports_lnk("report.md")+colors.RESET)


if __name__=="__main__":
    opt = ''

    tests_dir = os.path.realpath(os.path.dirname(__file__))#;print(tests_dir)
    hostname    = check_output('hostname -A', shell=1).decode().strip().split()[-1]
    reports_lnk = lambda file:"http://%s:8010/%s" %(hostname,file)
    run_opt = True
    libs = ['blochwave','EDutils']#,'multislice']
    if len(sys.argv)>1:
        opt = sys.argv[1]
        if 'h' in opt:
            print(usage )
            sys.exit()
        if 'f' in opt:
            reports_lnk=lambda file:"file://"+os.path.realpath(file)
        if 's' in opt:
            run_opt=False
        if 'l' in opt:
            libs=sys.argv[2].split(',')


    run_tests()
    # write_report()
