from mudata_explorer import app
from streamlit.delta_generator import DeltaGenerator
from plotly import graph_objects as go


def plot_mdata(
    container: DeltaGenerator,
    hspace=10,
    vspace=10
):
    if not app.has_mdata():
        return

    mdata = app.get_mdata()

    # Set up a figure
    fig = go.Figure()

    # Get the number of observations
    nobs = mdata.shape[0]

    # Plot the .obs as a filled rectangle with hovertext
    _add_box(
        fig,
        x=0,
        width=-mdata.obs.shape[1],
        y=0,
        height=-nobs,
        label="Observation Metadata",
        textangle=90
    )

    # Set up the current horizontal position for plotting
    hpos = hspace

    # For each modality
    for mod_name in mdata.mod.keys():

        # Get the number of features
        nvar = mdata.mod[mod_name].shape[1]

        # Plot the data (.X)
        _add_box(
            fig,
            x=hpos,
            width=nvar,
            y=0,
            height=-nobs,
            label=f"{mod_name} Features"
        )

        # Plot the metadata (.var)
        var_height = mdata.mod[mod_name].var.shape[1]
        _add_box(
            fig,
            x=hpos,
            width=nvar,
            y=vspace,
            height=var_height,
            label=f"{mod_name} Feature Annotations"
        )

        # Set up an offset for the .varm
        voffset = nobs + vspace

        # For each .varm
        for varm_key in mdata.mod[mod_name].varm.keys():

            # Get the number of columns
            ncol = mdata.mod[mod_name].varm[varm_key].shape[1]

            # Plot the .varm
            _add_box(
                fig,
                x=hpos,
                width=nvar,
                y=-voffset,
                height=-ncol,
                label=f"{mod_name} .varm[{varm_key}]"
            )

            # Update the offset
            voffset += vspace

        # For each .obsm
        hpos = hpos + nvar
        for obsm_key in mdata.mod[mod_name].obsm.keys():

            hpos += hspace

            # Get the number of columns
            ncol = mdata.mod[mod_name].obsm[obsm_key].shape[1]

            # Plot the .obsm
            _add_box(
                fig,
                x=hpos,
                width=ncol,
                y=0,
                height=-nobs,
                label=f"{mod_name} .obsm[{obsm_key}]",
                textangle=90
            )

            # Update the offset
            hpos += ncol

        # Increment the horizontal offset
        hpos += hspace

    # Turn off the axis ticks and labels
    fig.update_xaxes(showticklabels=False)
    fig.update_yaxes(showticklabels=False)

    # Remove the lines
    fig.update_layout(
        showlegend=False,
        hovermode="closest",
        plot_bgcolor="white",
        margin=dict(l=0, r=0, t=0, b=0),
        paper_bgcolor="white"
    )

    container.plotly_chart(fig)


def _add_box(
    fig: go.Figure,
    x: int,
    width: int,
    y: int,
    height: int,
    label: str,
    textangle=0
):

    if width == 0 or height == 0:
        return

    x1 = x + width
    y1 = y + height
    text = label + f"<br>({abs(height):,} rows x {abs(width):,} columns)"
    fig.add_trace(
        go.Scatter(
            x=[x, x1, x1, x, x],
            y=[y, y, y1, y1, y],
            mode="lines",
            fill="toself",
            text=text,
            hoverinfo="text",
            line=dict(width=0),
            fillcolor="rgba(173, 216, 230, 0.5)"  # Light blue color
        )
    )
    fig.add_annotation(
        x=x + (width / 2),
        y=y + (height / 2),
        text=text,
        textangle=textangle,
        showarrow=False
    )
