import altair as alt

# altair config
def sge_theme(*args, **kwargs):
    return {
        #'width': 800,
        #'height': 600,
        'config': {
            'title': {
                'align': 'center',
                'fontSize': 18,
                'labelFontSize': 15,
                'labelLimit': 0,
            },
            'axis': {
                'titleFontSize': 15,
                'labelFontSize': 15,
            },
            'legend': {
                'titleFontSize': 15,
                'labelFontSize': 16,
            }
        }
    }

# custom shape range gives us 10 distinct shapes
shape_range = [
    "circle",
    "square",
    "cross", 
    "diamond",
    "triangle-up",
    "triangle-down",
    "triangle-left",
    "triangle-right",
    "M0,.5L.6,.8L.5,.1L1,-.3L.3,-.4L0,-1L-.3,-.4L-1,-.3L-.5,.1L-.6,.8L0,.5Z",
    "arrow"
]

# enable altair theme
alt.data_transformers.enable('default', max_rows=None)
alt.themes.register('sge_theme', sge_theme)
alt.themes.enable('sge_theme')
