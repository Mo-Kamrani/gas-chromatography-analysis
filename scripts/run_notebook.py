"""Run this pure-Python notebook sequentially with headless Matplotlib.

This project-specific runner captures print output and plt.show() figures.
It does not emulate a Jupyter kernel or support IPython magics and rich displays.
"""

from contextlib import redirect_stderr, redirect_stdout
from io import BytesIO, StringIO
from pathlib import Path
import base64
import importlib.metadata
import os
import platform

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nbformat


def main():
    root = Path(__file__).resolve().parents[1]
    path = root / 'gas_chromatography.ipynb'
    notebook = nbformat.read(path, as_version=4)
    nbformat.validate(notebook)
    namespace = {'__name__': '__main__'}
    original_directory = Path.cwd()
    original_show = plt.show
    executed = 0
    try:
        os.chdir(root)
        for index, cell in enumerate(notebook.cells, start=1):
            if cell.cell_type != 'code':
                continue
            executed += 1
            cell.outputs = []
            cell.execution_count = executed
            stdout, stderr = StringIO(), StringIO()

            def flush_streams():
                for name, stream in (('stdout', stdout), ('stderr', stderr)):
                    text = stream.getvalue()
                    if text:
                        cell.outputs.append(nbformat.v4.new_output('stream', name=name, text=text))
                        stream.seek(0)
                        stream.truncate(0)

            def capture_figures(*args, **kwargs):
                flush_streams()
                for number in plt.get_fignums():
                    figure = plt.figure(number)
                    png = BytesIO()
                    figure.savefig(png, format='png', dpi=130, bbox_inches='tight')
                    cell.outputs.append(nbformat.v4.new_output(
                        'display_data', data={
                            'image/png': base64.b64encode(png.getvalue()).decode('ascii'),
                            'text/plain': f'Matplotlib figure {number}',
                        }, metadata={},
                    ))
                    plt.close(figure)

            plt.show = capture_figures
            with redirect_stdout(stdout), redirect_stderr(stderr):
                exec(compile(cell.source, f'{path.name}:cell-{index}', 'exec'), namespace)
                capture_figures()
            flush_streams()
    finally:
        plt.show = original_show
        plt.close('all')
        os.chdir(original_directory)
    notebook.metadata.language_info = {'name': 'python', 'version': platform.python_version()}
    notebook.metadata.execution_environment = {
        'engine': 'Sequential Python cells with headless Matplotlib (not a Jupyter kernel)',
        'python': platform.python_version(),
        **{name: importlib.metadata.version(name) for name in ('numpy', 'matplotlib', 'scipy', 'nbformat')},
    }
    nbformat.validate(notebook)
    temporary_path = path.with_suffix('.executed.ipynb')
    nbformat.write(notebook, temporary_path)
    temporary_path.replace(path)
    images = sum('image/png' in output.get('data', {}) for cell in notebook.cells for output in cell.get('outputs', []))
    print(f'Executed {executed} code cells; saved {images} embedded figures.')
    print(f'Notebook: {path.name}')
    print('PNG files: results/')


if __name__ == '__main__':
    main()
