function tout = plot_volcano_database(volcano,range,varargin)
nVarargin = length(varargin);
cd('~/observatorios/xlsx_databases/');
filename = [volcano,'_current.xlsx'];

figure;
hold on;
for i = 1:nVarargin
    sheet = varargin{i};
    disp(['Reading Sheet: ',sheet]);
    num = xlsread(filename,sheet);
    t = x2mdate(num(:,1)); %clobber old vectors (if they exist)
    t = datetime(t,'ConvertFrom','datenum');
    tI = t >= range(1) & t < range(2);
    t = t(tI);
    %t = datetime(t,'ConvertFrom','datenum');
    h(i) = plot(t,1:length(t));
    tout{i} = t;
end
legend(h,varargin,'location','northwest')
ax = gca;
ax.XTickLabelRotation = 45;
