%ang = angle(X(1,:)+1i*Y(1,:));
%r = sqrt(X(1,:).^2 + Y(1,:).^2);
%indBlue = find(r < 10 + 10/3);
%indGreen = find(r >= 10 + 10/3 & r < 10+20/3);
%indCyan = find(r >= 10+20/3);

theta = linspace(0,2*pi,1000);

ntime = size(X,1);
count=1;
%ntrac=4*floor(ntrac/4);
for i=1:1:ntime
  vecx = [posx(:,:,i) ;posx(1,:,i)];
  vecy = [posy(:,:,i) ;posy(1,:,i)];
  %plot(X(i,:),Y(i,:),'b.','markersize',1)
  %plot(X(i,:),Y(i,:),'b.')
  clf
  hold on
  if (options.confined)
    plot(10*exp(1i*theta),'k','linewidth',2);
    plot(20*exp(1i*theta),'k','linewidth',2);
  end

  %plot(X(i,indBlue),Y(i,indBlue),'b.')
  %plot(X(i,indGreen),Y(i,indGreen),'g.')
  %plot(X(i,indCyan),Y(i,indCyan),'c.')
  plot(X(i,:),Y(i,:),'b.')
  plot(vecx,vecy,'r-','linewidth',2)
  
  axis equal
  if (options.confined)
    axis([-21 21 -21 21]);
  else
    if (strcmp(file,'example1_data.bin'))
      axis([-3 3 -3 3])
    elseif (strcmp(file,'example2_data.bin'))
      axis([-3 3 -3 3])
    elseif (strcmp(file,'example3_data.bin'))
      axis([-3 12 -3 3])
    elseif (strcmp(file,'example4_data.bin'))
      axis([-3 3 -3 3])
    elseif (strcmp(file,'example6_data.bin'))
      axis([-0.1 pi+0.1 -0.1 pi+0.1])
    end
  end

  set(gca,'xtick',[])
  set(gca,'ytick',[])
  set(gca,'xcolor','w')
  set(gca,'ycolor','w')
  set(gca,'visible','on')
  %titleString = ['t = ' num2str(time(i),'%4.3f'),...
  %  ',\Delta A/A_0 = ' num2str(ea(i),'%4.3e'),...
  %  ',\Delta L/L_0 = ' num2str(el(i),'%4.3e')];
  titleString = ['t = ' num2str(time(i),'%4.3f')];
  title(titleString,'fontsize',16)
  filename = ['./image', sprintf('%04d',count), '.png'];
  count = count+1;
  figure(1);
  pause(0)
  %print(gcf,'-dpng','-r300',filename)
  hold off

end
