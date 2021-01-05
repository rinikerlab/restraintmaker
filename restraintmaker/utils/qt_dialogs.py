import os
import typing as t

import __main__
if not os.path.basename(__main__.__file__).startswith("test"):  #check, not in unit-tests
    from PyQt5 import QtWidgets

from restraintmaker.utils.Utilities import print


def create_input_dialog(parent_window, title: str) -> t.Callable[[str], str]:
    """
        create_input_dialog returns a function which:
    Opens a pyQt5 input Dialog with the speciefed message. In this format it can be provided to the get_args functions

    Parameters
    ----------
    parent_window
        parent QT window of the sub dialog

    title: str
        dialog title

    Returns
    -------
    t.Callable[str,[str])
        Returns a function, that opens an input dialog, with the specified title, and is centered at the interface_Pymol window.
    """

    def input_dialog_with_title(message: str)->t.Union[str, None]:
        """
            opens a dialog and returns the result as str

        Parameters
        ----------
        message: str
            Message to be dispayed in the input dialog

        Returns
        -------
        t.Union[str, None]
            input result

        """
        input, ok_pressed = QtWidgets.QInputDialog.getText(parent_window, title, message)  # self.cmd.gui.get_qtwindow()
        if ok_pressed:
            return input
        else:
            return None  # Do not catch the Error here. I'd prefer to do that in get args

    return input_dialog_with_title


def create_file_open_dialog():
    """
       create_file_dialog returns a function which:
       Opens a pyQt5 file Dialog with the speciefed message. In this format it can be provided to the get_args functions

    Returns
    -------
    t.Callable[str,[str])
        Returns a function, that opens an input dialog, with the specified title, and is centered at the interface_Pymol window.
    """

    def file_dialog_with_title(message: str):
        """
            opens a open-file dialog retrieving the correct path

        Parameters
        ----------
        message: str
            is ignored

        Returns
        -------
        str
            file path
        """
        path = QtWidgets.QFileDialog.getOpenFileName(None, 'Which File should be opened?', directory=os.getcwd())[0]
        return path

    return file_dialog_with_title


def create_file_save_dialog():
    """
       create_file_dialog returns a function which:
       Opens a pyQt5 file Dialog with the speciefed message. In this format it can be provided to the get_args functions

    Returns
    -------
     t.Callable[str,[str]]
        Returns a function, that opens an input dialog, with the specified title, and is centered at the interface_Pymol window.
    """

    def file_dialog_with_title(message: str):
        """
            opens a save-file dialog,getting the target path

        Parameters
        ----------
        message: str
            where to save?

        Returns
        -------
        str
            path to save to
        """
        path = QtWidgets.QFileDialog.getSaveFileName(None, message, directory=os.getcwd())[0]
        return path

    return file_dialog_with_title



def create_multi_dialog(title, inputs: t.List[str], options: t.Dict[str, t.List[str]] = {}, default: t.Dict = {}) -> \
t.Callable[[str], t.List[str or t.List[str]]]:
    '''
    create_multi_dialog return

    TODO get_args: I am not happpy with the create_multi_dialog solution: Now we have to 'prepare' the input function in set_x_type by checking what x(optimizer, filter...) type we want to create.
     The initial idea was that input_function is passed to get_args, which can then call it as needed, without having to know the format of its arguments
     We have to choose either:1) input function returns exactly one str. if more inputs are required, input fkt is called several times. => Does not allow to ask for all in one dialog
                       or    :2) Abolish get_args. Each type has a list required_args. the logic handler will read that list and then get them all from the user. The Optimizer?Filter... then only has to check the type

    :warning: Do not call before the constructor o f the Pymol module has finished, or pyumol will hate you.
    :warning: Will not throw an error if close button is pressed. Just return an enmpty list
    :param title: Title of theDialog
    :type title: str
    :param inputs; all inputs that should be asked for
    :type  inputs: t.List[str]
    :param options: optional, specifies the options an input can have. (input as key, list of acceptable values as values)
    :type options: t.dict[str,t.List[str]]
    :param default: Specifies the options an input can have. (input as key, list of acceptable values as value0
    :type default: t.dict[str,t.Any]
    :return: A function which will open a dialog, asking for the inputs specified in inputs
    :rtype: t.Callable[[str],t.List[str or t.List[str]]]
    '''

    class MultiDialog(QtWidgets.QDialog):
        'can be used to create input functions with as many arguments as necessary on the spot'

        def __init__(self):
            super().__init__()
            self.setWindowTitle(title)

        # def __del__(self):
        #     print('HELP! I am beeing deleted',mv=1)

    def _show_multi_dialog(dummy_str):
        '''The whole function _show_multi_dialog will be returned and can then be passed to get_args'''
        answers = []
        dial = MultiDialog()
        # In python components are added to the layot, not the window
        lay = QtWidgets.QGridLayout()

        # CREATE Components and add them to the layout
        fields: t.list[QtWidgets.QLineEdit or QtWidgets.QComboBox] = []
        for i_i, input in enumerate(inputs):
            if input in options.keys():
                new_field = QtWidgets.QComboBox()
                new_field.addItems(options[input])
                if input in default.keys():
                    new_field.setCurrentText(default[input])
            else:
                new_field = QtWidgets.QLineEdit()
                if input in default.keys():
                    new_field.setText(default[input])

            fields.append(new_field)  # TODO: Acutallz I can just loop over the grid layout to get all components
            lay.addWidget(QtWidgets.QLabel(input + ':'), i_i, 1)  # Label
            lay.addWidget(new_field, i_i, 2)

        OK_Button = QtWidgets.QPushButton('OK')

        def _ok_button_pressed():
            '''iterate over all input fields and add their value to answers'''
            for field in dial.children():
                # field = lay.itemAt(i)
                if isinstance(field, QtWidgets.QLineEdit):
                    answers.append(field.text())
                elif isinstance(field, QtWidgets.QComboBox):
                    answers.append(field.currentText())
                # else: QLabel etc => Do nothing
            dial.close()
            print(answers, mv=1)

        OK_Button.clicked.connect(_ok_button_pressed)

        lay.addWidget(OK_Button, len(inputs), 2)
        dial.setLayout(lay)
        dial.exec_()  # exec interrupts program flow unitl dialog is closed
        return answers

    return _show_multi_dialog