from bs4 import BeautifulSoup
import os,sys,argparse #.path import path
from subprocess import check_output,Popen,PIPE
from colorama import Fore as colors

usage='''Usage:
report.py [opt] [libs]

OPTIONS
-------
    opt  : h(help) s(skip running tests) l(libs) f(file serv) p(python)
    libs : <lib1>,<lib2>,.. if 'l' in opts
    python : python executable for subprocess
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
        cmd ='%s -m pytest %s' %(python,args)
        job ='cd %s/%s; if [ ! -d report ];then mkdir report;fi; %s ' %(tests_dir,lib,cmd)
        print(cmd)
        if run_opt:
            # oe=check_output(job, shell=1).decode().strip()
            p=Popen(job, shell=1);p.wait();o,e=p.communicate()
            if o  or e:print(o,e)

        #get summary
        lib_report=get_report_summary(report_file)
        lib_cov   =get_coverage_summary(cov_file)

        report+="%s | %s | %s \n" %(lib,lib_report,lib_cov)

    report_file=os.path.join(tests_dir,"report.md")
    with open(report_file,"w") as f: f.write(report)

    #convert to html
    html_file=report_file.replace('.md','.html')
    oe=check_output("pandoc %s > %s" %(report_file,html_file),shell=True).decode()
    print(colors.YELLOW+reports_lnk(html_file)+colors.RESET)


if __name__=="__main__":
    opt = ''

    tests_dir = os.path.realpath(os.path.dirname(__file__))#;print(tests_dir)
    hostname    = check_output('hostname -A', shell=1).decode().strip().split()[-1]
    reports_lnk = lambda file:"http://%s:8010/%s" %(hostname,file.split('tests/')[1])
    # libs = ['blochwave','EDutils']#,'multislice']
    # if len(sys.argv)>1:
    #     opt = sys.argv[1]
    #     if 'h' in opt:
    #         print(usage )
    #         sys.exit()
    #     if 'f' in opt:
    #         reports_lnk=lambda file:"file://"+os.path.realpath(file)
    #     if 's' in opt:
    #         run_opt=False
    #     if 'l' in opt:
    #         libs=sys.argv[2].split(',')
    #     if 'p' in opt:
    #         python = sys.argv[3]


    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--file'  ,action='store_true')
    parser.add_argument('-s','--skip'  ,action='store_true')
    parser.add_argument('-l','--libs'  ,default='blochwave,EDutils')
    parser.add_argument('-p','--python',default='python3')

    args = parser.parse_args()
    if args.file:reports_lnk=lambda file:"file://"+os.path.realpath(file)
    run_opt=not args.skip
    libs  = [ s.replace('/','') for s in args.libs.split(',')]
    python=args.python

    print(run_opt,libs,python)
    run_tests()
    # write_report()
