import pandas as pd


def aggregate_df_by_columns(df, key: str):
    """
    Aggregate pandas df along across rows/observations. 
    For numerical cols, do the mean, for categorical only keep those with homogeneous elements within the grouping key.

    Input:
        - df (pd.DataFrame)
        - key (str): column name according to which to group.

    Output:
        - grouped df
    """

    grouped = df.groupby(key)

    # Make dict: col -> agg func
    col_agg_func = {}
    
    for col in df.columns:

        if col == key:
            continue
        
        # Mean if numerical
        if pd.api.types.is_numeric_dtype(df[col]):
            
            col_agg_func[col] = 'mean'

        # If categorical, object or string ...
        elif pd.api.types.is_categorical_dtype(df[col]) or pd.api.types.is_object_dtype(df[col]) or pd.api.types.is_string_dtype(df[col]):
    
            # ... element if homogeneous, else None
            col_agg_func[col] = lambda x: x.iloc[0] if x.nunique() == 1 else None
    
        else:
            
            pass
            
    df_agg = grouped.aggregate(col_agg_func)


    # Exlude None categorical columns
    
    exclude = []
    
    for col in df_agg.columns:
    
        if df_agg[col].unique().all() == None:
    
            exclude.append(col)

    df_agg = df_agg.drop(columns=exclude)


    # Rename mean cols

    col_map = {}

    for col in df_agg.columns:

        if pd.api.types.is_numeric_dtype(df_agg[col]):

            col_map[col] = f'mean_{col}'

        else:
        
            col_map[col] = col

    df_agg = df_agg.rename(columns=col_map)


    return df_agg

