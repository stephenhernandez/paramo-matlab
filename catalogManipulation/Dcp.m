function D = Dcp(c,p,tStart,tStop)
if p == 1
    D = (1./c)./((log(1+(tStop./c))) - (log(1+(tStart./c))));
else
    D = ((1-p)./c)./(((1+(tStop./c)).^(1-p)) - ((1+(tStart./c)).^(1-p)));
end