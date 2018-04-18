function bold_fig
%
%     Make current plots publication bold  (Source unknown)
%     
%%orient landscape
orient portrait
hChildren0=get(gcf,'Children');
%% filter out uicontrol Children of matlab5
j=0;
for i=1:length(hChildren0)
    if strcmp(get(hChildren0(i),'Type'),'axes')
        j=j+1;
        hChildren(j)=hChildren0(i);
    end
end

%% back to matlab4 code:
nfigs=length(hChildren);

for i=1:nfigs
    hAxes=hChildren(i);
    
    set(hAxes,  'LineWidth',1.0);
    set(hAxes, 'FontWeight','bold');
    
    hXlabel =get(hAxes, 'xlabel');
    if exist('hXlabel')
        set(hXlabel,  'FontWeight', 'bold');
        %%set(hXlabel, 'FontSize', 14')
    end
    
    hYlabel = get(hAxes, 'ylabel');
    if exist('hYlabel')
        set(hYlabel,  'FontWeight', 'bold');
        %%set(hYlabel, 'FontSize', 14')
    end
    
    hZlabel = get(hAxes, 'zlabel');
    if exist('hZlabel')
        set(hZlabel,  'FontWeight', 'bold');
        %%set(hZlabel, 'FontSize', 14')
    end
    
    hTlabel = get(hAxes, 'title');
    if exist('hTlabel')
        set(hTlabel,  'FontWeight', 'bold')
        %%set(hTlabel, 'FontSize', 14')
    end
    
    hc = get(hAxes, 'Children');
    for i=1:length(hc)
        type = get(hc(i), 'type');
        if (length(type) ~=4)  %% skip over
        elseif (type=='text')
            set(hc(i), 'FontWeight', 'bold')
            %%set(hc(i), 'FontSize', 14)
        elseif (type=='line')
            set(hc(i), 'LineWidth', 2.0)
        end
    end
end
return
