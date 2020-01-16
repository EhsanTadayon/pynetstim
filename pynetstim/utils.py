### Utils
import os

def make_head_model(subject, freesurfer_dir):
    bem_folder = '{freesurfer_dir}/{subject}/bem'.format(freesurfer_dir=freesurfer_dir, subject=subject)
    os.makedirs(bem_folder,exist_ok=True)
    if not os.path.exists(os.path.join(bem_folder,'outer_skin_surface')):
        cmd ='cd {freesurfer_dir}/{subject}/bem; mri_watershed -surf surf {freesurfer_dir}/{subject}/mri/rawavg.mgz brain.mgz'.format(freesurfer_dir=freesurfer_dir, subject = subject)
        os.system(cmd)

        for f in ['lh.surf_brain_surface','lh.surf_inner_skull_surface','lh.surf_outer_skin_surface','lh.surf_outer_skull_surface']:
            cmd = 'mv {freesurfer_dir}/{subject}/bem/{f} {freesurfer_dir}/{subject}/bem/{f2}'.format(f=f,f2=f.split('lh.surf_')[1],subject=subject, freesurfer_dir=freesurfer_dir)
            os.system(cmd)
    else: 
        print('head model exists!')
        
        





### Python miscellaneous functions and classes for plotting 
### author: Ehsan Tadayon [ sunny.tadayon@gmail.com]
# Modified: June 17, 2019


def clean_plot(fig,ax,x_remove=['top'],y_remove=['right']):

    """ clean matplotlib plot
    
    parameters:
    ==========
    
    fig: matplotlib figure
    ax: matplotlib ax
    x_remove: the x axes that you want to remove ( either top or bottom)
    y_remove: the y axes that you want to remove (either rigth or left)
    
    
    returns:
    =======
    
    fig,ax 
    
    """
    
    x_keep = ['bottom','top']
    y_keep = ['left','right']
    
    
    
    for x_r in x_remove:
         x_keep.remove(x_r)
         
    for y_r in y_remove:
         y_keep.remove(y_r)
         
    try:
        for a in ax.flatten():
            
            for r in x_remove+y_remove:
                a.spines[r].set_visible(False)
            
            for k in x_keep:
                a.xaxis.set_ticks_position(k)
                
            for k in y_keep:
                a.yaxis.set_ticks_position(k)	
                
            a.tick_params(direction='out')
            
        return fig,ax
        
    except:
    
        for r in x_remove+y_remove:
            ax.spines[r].set_visible(False)
        
        for k in x_keep:
            ax.xaxis.set_ticks_position(k)
            
        for k in y_keep:
            ax.yaxis.set_ticks_position(k)
        
        ax.tick_params(direction='out')
        return fig,ax


# implementing classes to use to create simple html documents for quick report and visualization of images
# author: Ehsan Tadayon, M.D. [sunny.tadayo@gmail.com]
# modified: August 25, 2018

class HtmlDoc():

    def __init__(self,title,style=''):

        """
        html = HtmlDoc(title='test')
        html.add_header('Creating simple html document')
        html.add_image('test.png',500,500,'left')
        html.write('test.html')
        """
        
        self.html="""<!DOCTYPE html>
        <head>
            <title>{title}</title>
            <meta name="Ehsan Tadayon">
        <style>{style}</style>
            <body>
        """

        self.html=self.html.format(title=title,style=style)

    def add_style(self,style):

        pre=self.html.split('<style>')[0]
        post=self.html.split('<style>')[1]
        current_style=post.split('</style>')[0]
        post_style=post.split('</style>')[1]

        new_style=current_style+style
        self.html=pre+'<style>'+new_style+'</style>'+post_style


    def add_address(self,Author='Ehsan Tadayon',email='sunny.tadayon@gmail.com'):

        address="""
        <address>Author:{Author}<br>
        email:<a href="mailto:{email}">{email}</a></address>\n
        """

        address=address.format(Author=Author,email=email)
        self.html+=address


    def add_header(self,h,content):

        header='<{h}>{content}</{h}>\n'.format(h=h,content=content)
        self.html+=header


    def add_image(self,img_path,height,width,flt):

        img="""
        <img src="{img_path}" style="width:{width}px;height:{height}px;float:{flt};">\n"""
        img=img.format(img_path=img_path,height=height,width=width,flt=flt)
        self.html+=img

    def add_paragraph(self,paragraph):
        self.html+= paragraph+'\n'

    def add_tag(self,tag):
        self.html+=element.to_html()

    def write(self,fname):
        f=open(fname,'w')
        f.write(self.html)
        f.close()












            

    
 