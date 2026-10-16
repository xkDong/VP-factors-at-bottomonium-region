function [newsigma]=Newsigma(datacombis,l,curve)
m=length(datacombis(1,:));
newsigma=zeros(l,m);
sigma0=datacombis(2,:);
errsta=datacombis(3,:);
errsysunc=datacombis(4,:);
for i=1:l
    a=randn(1);b=randn(1);
    q=randn(1,m).*errsta(1,:);
    w=randn(1,m).*errsysunc(1,:);
    e=a*datacombis(5,:);
    newsigma(i,:)=sigma0+q+w+e+b*datacombis(6,:);
    for j=1:m
        while newsigma(i,j)<0
             newsigma(i,j)=sigma0(j)+randn(1).*errsta(j)+randn(1).*errsysunc(j)+randn(1)*datacombis(5,j)+randn(1)*datacombis(6,j);
        end
    end
end
newsigma(1,:)=sigma0;

% if curve==1
%     for i=1:l
%         newsigma(i,11:end)=smooth( datacombis(1,11:end),newsigma(i,11:end),0.03,'lowess');
%     end
% end

% %直线
% if curve==1
%     for i=1:m1
%         newsigma(:,i*k-k+1)=sigma1(:,i);
%     end
%     
%     for i=1:m1-1
%         %y=ax+b
%         dl=(newenergy(k*i-k+2)-newenergy(k*i-k+1));
%         a=(sigma1(:,i+1)-sigma1(:,i))./(newenergy(k*i+1)-newenergy(k*i-k+1));
%         for j=1:4
%             newsigma(:,k*i-k+1+j)=newsigma(:,k*i-k+j)+dl.*a;
%         end
%     end
% end

% %三次样条插值
% if curve==2
%     for i=1:l
%         sp=csape(energy,sigma1(i,:),'second',[0,0]);
%         newsigma(i,:)=fnval(sp,newenergy);
%         for j=1:m
%             if  newsigma(i,j)<0
%                 newsigma(i,j)=0.01;
%             end
%         end
%     end
% end
% 
% %平滑三次样条
% if curve==3
%     for i=1:l
%         sp=csaps(energy,sigma1(i,:));
%         newsigma(i,:)=fnval(sp,newenergy);
%         for j=1:m
%             if  newsigma(i,j)<0
%                 newsigma(i,j)=0.01;
%             end
%         end
%     end
% end
% 
% %B样条插值
% if curve==4
%     for i=1:l
%         sp=spapi(3,energy,sigma1(i,:));
%         newsigma(i,:)=fnval(sp,newenergy);
%         for j=1:m
%             if  newsigma(i,j)<0
%                 newsigma(i,j)=0.01;
%             end
%         end
%     end
% end
% 
% %均匀B样条
% if curve == 5
%     for i=1:l
%         values=spcrv([energy(1) energy energy(end);sigma1(i,1),sigma1(i,:),sigma1(i,end)],3);
%         L=length(values(1,:));
%         for j=1:m
%             for h=1:L
%                 if newenergy(j)==values(1,h)
%                     newsigma(i,j)=values(2,h);
%                     break;
%                 end
%                 if newenergy(j)<values(1,h)
%                     %直线插值
%                     newsigma(i,j)=(values(1,h)-values(1,h-1))*(newenergy(j)-values(1,h-1))/(values(1,h)-values(1,h-1))+values(2,h-1);
%                     break
%                 end
%             end
%         end
%         
%         for j=1:m
%             if  newsigma(i,j)<0
%                 newsigma(i,j)=0.01;
%             end
%         end
%     end
% end
end