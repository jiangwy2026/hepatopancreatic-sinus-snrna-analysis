import runpy, sys, os
log=open(r'D:\桌面\single cell\logs\sc_summary_new_names.log','w',encoding='utf8',buffering=1)
sys.stdout=log; sys.stderr=log
runpy.run_path(r'D:\桌面\single cell\analysis_output\scripts\176_sc_summary_new_names_preview.py', run_name='__main__')
log.close()
