from glob import glob
chans_dir = "../chans/"
stamps_dir = "../stamps/"
pv_dir = "../pv/"


f = open("translation_stamps.csv", 'w')

for i in range(42):
    nshell = i + 1
    if nshell < 10:    
        fig_str = "f4_0"+str(nshell)+".pdf"
    else:
        fig_str = "f4_"+str(nshell)+".pdf"
    
    title_str = "Shell " + str(nshell)
    f_cap = open(stamps_dir + fig_str[:-3] + "txt", 'r')
    caption_str = f_cap.read()
    print(caption_str)
    
    trans_str = '{}, {}, "{}"\n'.format(fig_str,title_str,caption_str)
    f.write(trans_str)
f.close()