function [Population] = FindParticles(Image,threshold,ROI_size_gauss,folder)
% feed me an image with AuNPs!!!!

if isempty(threshold) == 1
    threshold = [ 0.005 0.75 ]; % The thresholds are 0.005 and 0.75
else
end    
 
% ROI_size_gauss =ROI_size_gauss;

% Remove background
[Corrected, Background] = BackgroundCorrection_WaveletSet(Image,4);
Image = Image - Background;

% Find particles
load('Typical_BeadFilter_for_Finding_AuNPs.mat');
Population = LocateBeadsProfile(Image, BeadFilter, false, threshold(1), threshold(2),20);
Objects = Population.Location;

      figure
      imagesc(Image,[0 1e4])
      colorbar; colormap jet;
      for i_Beads = 1:Population.N_Object
          rectangle('Position',[Objects(i_Beads,1)-0.5*(ROI_size_gauss-1),Objects(i_Beads,2)-0.5*(ROI_size_gauss-1),ROI_size_gauss,ROI_size_gauss],'EdgeColor','white')
          text(Objects(i_Beads,1)+0.5*(ROI_size_gauss+3),Objects(i_Beads,2)-1,num2str(i_Beads),'Color','white','FontSize',8,'FontWeight','bold');
      end   
      title({'The localized particles'});
      print(gcf,[folder 'The localized particles.png'],'-dpng','-r300','-opengl') %save file as png
      saveas(gcf,[folder 'The localized particles.fig'],'fig') %save file as matlab fig




