import numpy as np
import pandas as pd


def get_overlapping_interval_indices(start, stop, start_array, stop_array, handle_overlapping):
    """Return the indices of all intervals stored in start_array, stop_array that overlap start, stop"""
    # TODO: replace with better algorithm - e.g. rustbio intervl trees
    if not handle_overlapping:
        first_start_smaller = np.searchsorted(start_array, start) - 1
        first_end_larger = np.searchsorted(stop_array, stop, "right") + 1
    else:
        first_start_smaller = 0
        first_end_larger = len(stop_array)
    result = []
    for possible_match in range(
        max(0, first_start_smaller), min(first_end_larger, len(stop_array))
    ):
        s = max(start, start_array[possible_match])
        e = min(stop, stop_array[possible_match])
        if s < e:
            result.append(possible_match)
    return result


def merge_intervals(df):
    """take a {chr, start, end, *} dataframe and merge overlapping intervals"""
    if hasattr(df, "to_pandas"):
        raise ValueError("pydataframe passed")
    df = df.sort_values(["chr", "start"], ascending=[True, True]).reset_index(
        drop=True
    )  # you need to do this here so it's true later...
    keep = _merge_choose_rows_to_keep_and_update_positon(df)
    return df.loc[keep].reset_index(drop=True)


def merge_intervals_with_callback(df, callback):
    """take a {chr, start, end, *} dataframe and merge overlapping intervals, calling callback for group larger than one.."""
    df = df.sort_values(["chr", "start"], ascending=[True, True]).reset_index(
        drop=True
    )  # you need to do this here so it's true later...
    keep = _merge_choose_rows_to_keep_and_update_positon(df)
    # now, I have each run identified by ending in keep=True, so it's easy to just walk them and call the merge function appropriatly
    # this saves me from having to figure out where the run started, which the normal merge doesn't need...
    last_keep = -1
    res = []
    for ii, do_keep in enumerate(keep):
        if do_keep:  # the end of a run...
            if last_keep < ii - 1:  # more than one...
                subset = df.iloc[last_keep + 1 : ii + 1]
                row_data = callback(subset)
                if not isinstance(
                    row_data, dict
                ):  # and not (isinstance(row_data, pd.core.series.Series) and len(row_data.shape) == 1):
                    print("type", type(row_data))
                    # print 'len(shape)', len(row_data.shape)
                    print(callback)
                    raise ValueError(
                        "Merge_function returned something other than dict (writing to the pandas series directly is very slow, call to_dict() on it, then modify it.)"
                    )
                if set(row_data.keys()) != set(df.columns):
                    raise ValueError(
                        "Merge_function return wrong columns. Expected %s, was %s"
                        % (df.columns, list(row_data.keys()))
                    )
                row_data["start"] = df.at[ii, "start"]
                row_data["stop"] = df.at[ii, "stop"]
                res.append(row_data)
            else:
                res.append(df.iloc[last_keep + 1].to_dict())
            last_keep = ii
    print(res)
    res = pd.DataFrame(res)[df.columns].reset_index(drop=True)
    return res


def _merge_choose_rows_to_keep_and_update_positon(df):
    """A helper for the merge_intervals and merge_intervals_with_callback functions that refactors some common code.
    Basically, updates continuous 'runs' of overlapping intervals so that the last one has start = start_of_first, end = end of longest,
    and then returns a boolean array with True if it's the last one"""
    if hasattr(df, "to_pandas"):
        raise ValueError("pydataframe passed")
    chrs = np.array(df["chr"])
    starts = np.array(df["start"])
    stops = np.array(df["stop"])
    ii = 0
    lendf = len(df)
    keep = np.zeros((lendf,), dtype=np.bool)
    last_chr = None
    last_stop = 0
    last_row = None
    while ii < lendf:
        if chrs[ii] != last_chr:
            last_chr = chrs[ii]
            last_stop = 0
        if starts[ii] < last_stop:
            starts[ii] = starts[last_row]
            stops[ii] = max(stops[ii], last_stop)
        else:
            if last_row is not None:
                keep[last_row] = True
                # new_rows.append(df.get_row(last_row))
        if stops[ii] > last_stop:
            last_stop = stops[ii]
        last_row = ii
        ii += 1
    if last_row is not None:
        keep[last_row] = True
    df["start"] = starts
    df["stop"] = stops
    return keep


def merge_identical(sub_df):
    if (len(sub_df["start"].unique()) != 1) or (len(sub_df["stop"].unique()) != 1):
        raise ValueError(
            "Overlapping intervals: %s, merge_identical was set." % (sub_df,)
        )
    return sub_df.iloc[0].to_dict()
