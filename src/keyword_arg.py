def set_defaults(kwargs_dict,defaults_dict):
    for pair in defaults_dict.items():
        kwargs_dict.setdefault(*pair)
    return kwargs_dict

def parse_keywords(arg_dict,prefixes):

    # Initialize output
    new_dicts = {}
    output_types = list(prefixes) + ['other']
    for type in output_types:
        new_dicts[type] = {}
        
    # Need to sort prefixes list before running the match
    # This ensures that arguments like 'psth_plot_kw'
    # get parsed *before* 'psth_kw'; otherwise, output would include
    # 'plot_kw' in the 'psth' dict instead of 'kw' in the 'psth_plot' dict.
    prefixes.sort(reverse=True)
    
    # Assign input arguments to separate dicts based on prefixes
    for arg in arg_dict:
        match = False
        for prefix in prefixes:
            if arg.startswith(prefix):
                new_kw = arg.replace(prefix,'',1).lstrip('_')
                new_dicts[prefix][new_kw] = arg_dict[arg]
                match = True
        if not match:
            new_dicts['other'][arg] = arg_dict[arg]
        
    return [new_dicts[type] for type in output_types]