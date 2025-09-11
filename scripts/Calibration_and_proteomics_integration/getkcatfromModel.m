function kcat = getkcatfromModel(model, reaction_list)
for i = 1:numel(reaction_list)
    reaction_name = reaction_list{i};
    reaction_index = find(strcmp(reaction_name, model.ec.rxns));
    kcat(i) = model.ec.kcat(reaction_index);
end
kcat = kcat';