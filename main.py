# # #  Global constants
main_screen_size = (750, 800)

from kivy.config import Config
Config.set('graphics', 'width', main_screen_size[0])
Config.set('graphics', 'height', main_screen_size[1])
Config.set('graphics', 'resizable', True)
Config.set('input', 'mouse', 'mouse,disable_multitouch')

# All config categories initialized in AlnKeyApp.build_config()
config_paths_category = 'User paths'
config_class_colours_category = 'Residue class colours'
config_variation_colours_category = 'Variation screen colours'
config_category_keys = [config_paths_category, config_class_colours_category, config_variation_colours_category]

import time
import threading

import os, sys
from app_scripts import sequ, align
from kivy.resources import resource_add_path
from kivy.app import App
from kivy.lang import Builder
from kivy.core.window import Window
from kivy.properties import StringProperty, BooleanProperty
from kivy.uix.screenmanager import ScreenManager, Screen
from kivy.uix.relativelayout import RelativeLayout
from kivy.uix.gridlayout import GridLayout
from kivy.uix.textinput import TextInput
from kivy.uix.dropdown import DropDown
from kivy.uix.button import Button
from kivy.uix.popup import Popup
from kivy.uix.label import Label
from kivy.uix.image import Image
from kivy.animation import Animation

from kivy.clock import Clock

from kivy.graphics import Color, Rectangle, Line
from kivy.core.text import Label as CoreLabel
from kivy.graphics.instructions import Callback


# BUGS
# draw_loading_graphic(self, dt=None) is messed up. It works great if you load a large aln. But only the first load. The graphic isn't drawn on subsequent loads. Unless you comment out the draw_graphics() call; somehow this makes it all work again. Right now it's using Clock.schedule_once, but I'd rather not. Once fixed, remove the 'dt=None' args from various functions.

# TODO
# Implement a thread pool for the App. Will probably fix the loading graphic bug.
# I want some methods to label regions. One line of graphics that is very general, for lobes, active sites, whatever. Another line specifically for secondary structure. Users can manually input ss boundaries (can then save as dssp file), or load one of the file formats used by DSSP (https://swift.cmbi.umcn.nl/gv/dssp/index.html); that page also shows how to use their API, so I should implement that (user loads PDB structure, sends to web server, waits for response).
#  - I like the graphics used in SSDraw (https://pmc.ncbi.nlm.nih.gov/articles/PMC10680343/)
# Reformat the landing page. I want one area (info/button) for an alignment file, a second for a PDB file, a third for secondary structure, fourth for regions, fifth for sequence groups.
#  - Depending on the files currently loaded, the other screens will become available. Ex only aln loaded, can view, graph quality, but can't map to pdb. If only pdb loaded, can only predict secondary structure. If only pdb + dssp loaded, can View Alignment (where the alignment is just the sequence from the pdb).
# Give the variation graph a y-axis title, possibly x-axis as well
# Consider sliding window for calculation of variation. Or perhaps calculation of the mean line. Check a couple real gappy alignments. 
# Implement ability to select one sequence to display instead of consensus. Searchable + scrollable popup list. The View screen will likely use a version that can incorporate checkboxes. 
#  - When displaying sequence on Variation, ignore gaps? Option to? 

# not sure if my .ico has multiple sizes; either way look into packaging different images for the different sizes into one .ico
# Probably do some optimization of the hover message framework. I can probably make the lookup much faster, there must be some data structure that's well suited to checking if a value falls within a selection of ranges.
#  - This is for get_hover_message()

# Could be cool to have a function to generate a sequence motif logo. Let the user select the boundaries, should be easy enough to do.

# I think I want a second draw option for alignments, essentially Figure 2 from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4155614/ looks like it would be a good way to see subtypes.
#  - Could even make that one interactive; click a line to get the names of all sequences making that transition.
#  - I really, really like the idea of clustering the Sequence Diversity Diagrams. Just k-medoids(means?), clustering the abstractions (line graphs) themselves. This will point out subtypes in the alignment. Or maybe EM? No reason to expect group sizes to be the same size.
# - Ensure there's a "gap" category.
# - Existing tool: https://www.ebi.ac.uk/research/goldman/software/alvis/
# - Related: this tool has a network diagram of sequence alignments which I like (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0259726)

# build options menu as described https://kivy.org/doc/stable/api-kivy.app.html

# For publishing:
# - My diversity diagram is described in the 2014 paper (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4155614/) but is actually a figure generated by "Parallel sets" (https://pubmed.ncbi.nlm.nih.gov/16805264/). The 2014 paper's work is not good, I'm not using it (it's also no longer available). I'm implementing a sequence alignment-specific adaptation of ParSet.
# - I probably need to evaluate/compare my variation score in some manner. This paper suggests ways to do that https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-7-484


# # #  Common accessory methods  # # #
def pick_nice_interval(number):
    """Picks a "nice" value that goes into the given number 4-10 times, where the value is a multiple of 10 times: 10, 20, or 50. The lowest returned value is 10. Made to help with labeling graph axes."""
    low_num, high_num = number//10, number//4
    mult_exp = max(len(str(low_num)) - 2, 0)
    mult = 10**mult_exp
    vals = {10}
    for val in (10, 20, 50):
        prod1, prod2 = val*mult, val*mult*10
        if low_num <= prod1 <= high_num:
            vals.add(prod1)
        if low_num <= prod2 <= high_num:
            vals.add(prod2)
    return max(vals)

# # #  Common elements  # # #
class FileChooserPopup(Popup):
    # Descendants must specify a `selection_made()` function.
    # On Windows, the console throws a ton of apparently harmless errors because building the file tree involves trying to access protected files. This is solved with "filters: [lambda dir, fname: not (fname.endswith('.sys') or fname.endswith('.tmp'))]" under the FileChooserListView element in the .kv
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.prev_path = ''
        # #  Added in aln_key_layout.kv
        # self.filechooser = FileChooserListView
        self.filechooser.bind(path=self.path_changed)
    def path_changed(self, filechooser, path):
        # Only saved if the user entered a real path
        self.prev_path = path
    def validate_user_path(self, textinput):
        new_path = textinput.text
        if os.path.isdir(new_path):
            self.filechooser.path = new_path
        elif os.path.isfile(new_path):
            self.filechooser.path = os.path.dirname(new_path)
            self.filechooser.selection = [new_path]
        else:
            textinput.text = self.prev_path
    def up_level(self):
        prev, _ = os.path.split(self.filechooser.path)
        self.filechooser.path = prev
    def open(self, title, init_path, callback):
        """The callback function will be called as `callback(dir_path, filepath/filename)` when the user presses "Load" or "Save". The callback will not be called if the user presses "Cancel"."""
        try:
            self.filechooser._update_files() # May not be stable, not officially supported. Used to force an update for new files. Hence the "try".
        except:
            pass
        self.title = title
        self.filechooser.path = init_path
        self.callback = callback
        super().open()
class LoadFilePopup(FileChooserPopup):
    # #  Attributes added in aln_key_layout.kv
    # self.loadbutton = Button
    # self.selectedlabel = Label
    def selection_made(self, selection):
        if selection:
            self.loadbutton.disabled = False
            path = selection[0]
            filename = os.path.basename(path)
        else:
            self.loadbutton.disabled = True
            filename = ''
        self.selectedlabel.text = filename
class SaveFilePopup(FileChooserPopup):
    # #  Attributes added in aln_key_layout.kv
    # self.savefilename = TextInput
    def selection_made(self, selection):
        if selection:
            path = selection[0]
            filename = os.path.basename(path)
            self.savefilename.text = filename

class PermissionPopup(Popup):
    message_text = StringProperty('')
    def __init__(self):
        super().__init__()
        self.ok_callback = None
        self.ok_args = ()
        self.cancel_callback = None
        self.cancel_args = ()
    def open(self, title, message, ok_callback=None, ok_args=(), cancel_callback=None, cancel_args=()):
        self.title = title
        self.message_text = message
        self.ok_callback, self.ok_args = ok_callback, ok_args
        self.cancel_callback, self.cancel_args = cancel_callback, cancel_args
        super().open()
    def ok_button_callback(self):
        if self.ok_callback:
            self.ok_callback(*self.ok_args)
    def cancel_button_callback(self):
        if self.cancel_callback:
            self.cancel_callback(*self.cancel_args)

class HeaderLogo(Image):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.source = os.path.join('app_data', 'icon.png')

class BaseDropDown(DropDown):
    def __init__(self, parent_screen, **kwargs):
        super().__init__(**kwargs)
        self.screen = parent_screen
        self.container.spacing = self.spacing # Normally can't really
        self.container.padding = self.padding # set these settings.
class AlnExportDropDown(BaseDropDown):
    pass
class AlnSelectChoice(BaseDropDown):
    """A DropDown object allowing the user to select one of a number of choices, storing the string value to self.choice. Use the add() function to add a choice to the list, and the clear() function to remove all current choices and reset self.choice back to the default_choice. You can pass a function as on_change, which will be called as fxn(new_choice) when the user selects a choice."""
    # I want to get a scrollbar working, maybe using a ScrollView?
    # Should probably try and replace DropDown.container with a RecycleView, as it can handle large number of items.
    def __init__(self, *args, on_change=None, default_choice='', **kwargs):
        self.choice = default_choice
        self.default_choice = default_choice
        if on_change:
            self.on_change = on_change
        super().__init__(*args, **kwargs)
    def add(self, text):
        # Should deal with updating the max width if 'text' is longer than the current max string
        self.add_widget(SelectChoiceEntry(text=text))
    def clear(self):
        self.clear_widgets()
        self.choice = self.default_choice
    def on_change(self, new_choice):
        pass
class SelectChoiceEntry(Button):
    pass

class ObjectGroupLayout(GridLayout):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    def do_layout(self, *args, **kwargs):
        # Dynamic heights were being a nightmare without this
        children_height = sum(child.height+self.spacing[1] for child in self.children if child.height>0) - self.spacing[1]
        if children_height > 0:
            self.height = children_height + self.padding[1] + self.padding[3]
        else:
            self.height = 0
        super().do_layout(*args, **kwargs)

class LabelledInput(RelativeLayout):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        #self.text_input = ExpandingTextInput
    def reset(self):
        self.text_input.reset()
    @property
    def text(self):
        return self.text_input.text
    @text.setter
    def text(self, val):
        self.text_input.text = val
class LabelledInputRange(RelativeLayout):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        #self.text_input1 = ExpandingTextInput
        #self.text_input2 = ExpandingTextInput
    def reset(self):
        self.text_input1.reset()
        self.text_input2.reset()
    @property
    def texts(self):
        return [self.text_input1.text, self.text_input2.text]
    @texts.setter
    def texts(self, vals):
        self.text_input1.text = vals[0]
        self.text_input2.text = vals[1]
    @property
    def text1(self):
        return self.text_input1.text
    @text1.setter
    def text1(self, val):
        self.text_input1.text = val
    @property
    def text2(self):
        return self.text_input2.text
    @text2.setter
    def text2(self, val):
        self.text_input2.text = val
class ExpandingTextInput(TextInput):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    def reset(self):
        self.text = self.default_input_text
    def on_focus(self, instance, has_focus):
        duration = self.parent.expand_duration
        if has_focus:
            ani = Animation(width=self.parent.expanded_width, duration=duration)
            ani.start(self)
        else:
            self._validate_format_text()
            self.cursor = (0, 0)
            ani = Animation(width=self.parent.input_width, duration=duration)
            ani.start(self)
    # # #  Private methods  # # #
    def _validate_format_text(self):
        fmt = self.parent.input_format
        text = self.text.strip()
        if text == '' and not self.parent.allow_empty:
            return self.reset()
        elif fmt == 'int':
            try:
                val = int(text)
            except:
                return self.reset()
            if self.min_value != None and val < self.min_value:
                return self.reset()
            if self.max_value != None and val > self.max_value:
                return self.reset()
            text = str(val)
        elif fmt == 'float':
            try:
                val = float(text)
            except:
                return self.reset()
            if self.min_value != None and val < self.min_value:
                return self.reset()
            if self.max_value != None and val > self.max_value:
                return self.reset()
        elif fmt == 'upper':
            text = text.upper()
        elif fmt == 'lower':
            text = text.lower()
        self.text = text


# # #  Inheritable base classes  # # #
class BaseScreen(Screen):
    """Should be inherited by all screen classes. Provides a few basic attributes."""
    def __init__(self, *args, **kwargs):
        self.app = App.get_running_app()
        self.screen_size = main_screen_size
        super().__init__(*args, **kwargs)
    def alignment_loaded(self):
        """Called by the MainScreen whenever a new alignment has been successfully loaded. Overwrite this to implement any sort of desired cache behaviours."""
        pass
class LoadFiles:
    """Provides the functionality to handle loading files, and managing the saving of those paths between program sessions.
    Expects that children define:
    self.path_config_keys - A list of strings, each config key for saving & loading files. Must be defined before the child calls `super().__init__()`

    Provides the following attributes:
    self.get_load_filepath(title, config_key, callback, callback_args=()) - Calls `callback(args)` where `args=(dirpath, filepath)+callback_args` once the user selects a file.
    """
    def __init__(self, *args, **kwargs):
        for key in self.path_config_keys:
            self.app.config_path_keys.add(key)
        self._load_config_key = None
        self._load_callback = None
        self._load_callback_args = ()
        super().__init__(*args, **kwargs)
    def get_load_filepath(self, title, config_key, callback, callback_args=()):
        self._load_config_key = config_key
        self._load_callback = callback
        self._load_callback_args = callback_args
        default_path = self.app.config.get(config_paths_category, config_key)
        self.manager.load_popup.open(title, default_path, self._load_file_chosen)
    def _load_file_chosen(self, dirpath, filepath):
        self.app.config.set(config_paths_category, self._load_config_key, dirpath)
        args = (dirpath, filepath) + self._load_callback_args
        self._load_callback(*args)
class SaveFiles:
    """Provides the functionality to handle saving files, and managing the saving of those paths between program sessions.
    Expects that children define:
    self.path_config_keys - A list of strings, each config key for saving & loading files. Must be defined before the child calls `super().__init__()`

    Provides the following attributes:
    self.get_save_filepath()
    """
    def __init__(self, *args, **kwargs):
        for key in self.path_config_keys:
            self.app.config_path_keys.add(key)
        self._save_suffix = None
        self._save_config_key = None
        self._save_callback = None
        self._save_callback_args = ()
        super().__init__(*args, **kwargs)
    def get_save_filepath(self, title, default_name, suffix, config_key, callback, callback_args=()):
        self._save_suffix = suffix
        self._save_config_key = config_key
        self._save_callback = callback
        self._save_callback_args = callback_args
        default_path = self.app.config.get(config_paths_category, config_key)
        self.manager.save_popup.savefilename.text = default_name
        self.manager.save_popup.open(title, default_path, self._save_file_chosen)
    def _save_file_chosen(self, dirpath, filename):
        if filename:
            if self._save_suffix and not filename.lower().endswith(self._save_suffix):
                filename += self._save_suffix
            if os.path.isfile(os.path.join(dirpath, filename)):
                perm_msg = 'The destination:\n"{}"\n\nalready has a file named:\n"{}"\n\nOverwrite and replace it?'.format(dirpath, filename)
                self.manager.permission_popup.open('File overwrite warning', perm_msg, ok_callback=self._do_save_callback, ok_args=(dirpath, filename))
            else:
                self._do_save_callback(dirpath, filename)
    def _do_save_callback(self, dirpath, filename):
        self.app.config.set(config_paths_category, self._save_config_key, dirpath)
        args = (dirpath, filename) + self._save_callback_args
        self._save_callback(*args)
class DrawGraphics:
    """Base class inherited by X, Y, Z. Provides methods to draw on a given canvas.

    Children are expected to define:
    self.draw_canvas = BoxLayout
    self.draw_canvas_view = ScrollView

    Optionally children can define/overwrite:
    self.screen_size = (width, height)
    """
    def __init__(self, *args, **kwargs):
        #super().__init__(*args, **kwargs)
        self.registered_colour_dicts = {}
        self.residue_classes = {
            'E':'acidic', 'D':'acidic',
            'N':'polar', 'Q':'polar',
            'R':'basic', 'K':'basic', 'H':'basic',
            'C':'cysteine',
            'M':'hydrophobic', 'V':'hydrophobic', 'L':'hydrophobic', 'I':'hydrophobic',
            'Y':'large', 'F':'large', 'W':'large',
            'A':'small', 'S':'small', 'T':'small',
            'P':'turn', 'G':'turn'}
        default_class_colours = {'acidic':'255,0,0,255', 'polar':'161,255,10,255', 'basic':'20,125,245,255', 'cysteine':'222,255,10,255', 'hydrophobic':'255,211,0,255', 'large':'255,135,0,255', 'small':'10,255,153,255', 'turn':'190,10,255,255'}
        # unused: cyan (10/255, 239/255, 255/255, 1), dark purple (88/255,10/255,255/255,1)
        self.register_colour_dict(config_class_colours_category, default_class_colours)
        self.class_colours = self.get_colour_dict(config_class_colours_category)
    # #  Common methods
    def get_residue_colour(self, residue):
        """Returns a tuple of floats (int / 255), or None if the residue is not found."""
        residue_class = self.residue_classes.get(residue, None)
        return self.class_colours.get(residue_class, None)
    def register_colour_dict(self, config_key, default_colours):
        """Colours in the default dict should be of the form {'key':'r,g,b,a', ...} where r/g/b/a are values from 0-255. Loads the colours from the config file under the given key, only resorting to the default if the key is not found or is malformed."""
        cols = self.load_config_colour_dict(config_key, default_colours)
        self.registered_colour_dicts[config_key] = {'dict':cols, 'defaults':default_colours}
    def get_colour_dict(self, config_key):
        """Returns a registered colour dict. Values in that dict should not be manually changed, instead use self.set_colour_dict(config_key, key, r,g,b,a), or self.reset_colour_dict(config_key, key)."""
        return self.registered_colour_dicts[config_key]['dict']
    def set_colour_dict(self, config_key, key, r,g,b,a):
        """r/g/b/a are values from 0-255."""
        for i in (r,g,b,a):
            if not 0<=i<=255:
                raise ValueError('error setting colour dict. The given rgba values must be between 0 and 255.') 
        cols = self.get_colour_dict(config_key)
        colour_tup = (r/255, g/255, b/255, a/255)
        cols[key] = colour_tup
        colour_str = '{},{},{},{}'.format(r,g,b,a)
        self.app.config.set(config_key, key, colour_str)
    def reset_colour_dict(self, config_key, key=None):
        """Resets the key for a registered colour dict to the original default. If key=None, resets all keys."""
        cols = self.get_colour_dict(config_key)
        defs = self.registered_colour_dicts[config_key]['defaults']
        if not key:
            keys = defs.keys()
        else:
            keys = [key]
        for k in keys:
            r,g,b,a = map(float, defs[k].split(','))
            self.set_colour_dict(config_key, k, r,g,b,a)
    # # #  Drawing methods
    def draw_methods_factory(self):
        # TODO: rename to get_draw_methods
        """Returns the drawing methods as a dict: {'draw_rect':fxn, 'draw_label':fxn}. Allows the returned drawing methods to use purely relative coordinates. The defined methods:
        draw_line()
        draw_curve()
        draw_rect()
        draw_label()
        """
        pos_x, pos_y = self.draw_canvas.pos
        def draw_line(rel_points, width=1.0, colour=None, dashes=False):
            # rel_points should be a list of tuples: [(x1,y1), (x2,y2), ...] in relative coordinates
            # dashes should be (dash_len, space_len) to be dashed; also requires that width=1.0
            points = [(pos_x+x, pos_y+y) for x,y in rel_points]
            if colour:
                Color(rgba=colour)
            if dashes != False and width == 1.0:
                dash_len, space_len = dashes
                return Line(points=points, width=width, cap='none', dash_length=dash_len, dash_offset=space_len)
            else:
                return Line(points=points, width=width, cap='none')
        def draw_curve(rel_points, width=1.1, colour=None, continue_from=False, continue_to=False):
            """Generates 2 new points between each point in rel_points to make a smoother curve. If continue_from / continue_to are given as points, they won't be drawn but the curve at the start / end will bend towards them."""
            def midpoints(prev_pnt, cur_pnt, cont_pnt=0):
                """Generates 2 points between prev_pnt and cur_pnt. These points make an approximation of a curve because dx and dy are applied in different proportions. cont_pnt is used to indicate if midpoint 1 or 2 is towards a continuation point."""
                d_x, d_y = cur_pnt[0] - prev_pnt[0], cur_pnt[1] - prev_pnt[1]
                if cont_pnt == 1:
                    mid1 = (prev_pnt[0] + d_x/4, prev_pnt[1] + d_y/2)
                else:
                    mid1 = (prev_pnt[0] + d_x/4, prev_pnt[1] + d_y/6)
                if cont_pnt == 2:
                    mid2 = (prev_pnt[0] + d_x*3/4, prev_pnt[1] + d_y/2)
                else:
                    mid2 = (prev_pnt[0] + d_x*3/4, prev_pnt[1] + d_y*5/6)
                return mid1, mid2
            points = [(rel_points[0][0]+pos_x, rel_points[0][1]+pos_y)]
            prev_pnt = points[0]
            for x, y in rel_points[1:]:
                cur_pnt = (x+pos_x, y+pos_y)
                mid1, mid2 = midpoints(prev_pnt, cur_pnt)
                points.append(mid1)
                points.append(mid2)
                points.append(cur_pnt)
                prev_pnt = cur_pnt
            if continue_from:
                mid1, mid2 = midpoints((continue_from[0]+pos_x, continue_from[1]+pos_y), points[0], cont_pnt=1)
                points.insert(0, mid1)
                points.insert(1, mid2)
            if continue_to:
                mid1, mid2 = midpoints(points[-1], (continue_to[0]+pos_x, continue_to[1]+pos_y), cont_pnt=2)
                points.append(mid1)
                points.append(mid2)
            if colour:
                Color(rgba=colour)
            return Line(points=points, width=width, cap='round', joint='bevel')
        def draw_rect(pos, size, box_colour=None):
            if box_colour:
                Color(rgba=box_colour)
            return Rectangle(pos=(pos_x+pos[0], pos_y+pos[1]), size=size)
        def draw_label(pos, text, size=None, halign='left', valign='bottom', padding_x=0, padding_y=0, font_size=12, font_colour=None, box_colour=None):
            lab = CoreLabel(text=text, font_size=font_size, bold=True)
            lab.refresh()
            lab_w, lab_h = lab.texture.size
            if size == None:
                size = (0, 0)
            box_size = max(lab_w+padding_x, size[0]), max(lab_h+padding_y, size[1])
            box_pos = round(pos_x+pos[0]), round(pos_y+pos[1])
            text_pos = round(box_pos[0]+box_size[0]/2.0-lab_w/2.0), round(box_pos[1]+box_size[1]/2.0-lab_h/2.0)
            x_off, y_off = 0, 0
            if halign == 'middle':
                x_off = box_size[0] / 2
            elif halign == 'right':
                x_off = box_size[0]
            if valign == 'middle':
                y_off = box_size[1] / 2
            elif valign == 'top':
                y_off = box_size[1]
            box_pos = (box_pos[0]-x_off, box_pos[1]-y_off)
            text_pos = (text_pos[0]-x_off, text_pos[1]-y_off)
            if box_colour:
                Color(rgba=box_colour)
                Rectangle(pos=box_pos, size=box_size)
            if font_colour:
                Color(rgba=font_colour)
            rec = Rectangle(texture=lab.texture, pos=text_pos, size=(lab_w, lab_h))
            del lab
            return rec
        return {'draw_line':draw_line, 'draw_curve':draw_curve, 'draw_rect':draw_rect, 'draw_label':draw_label}
    # # #  Private methods
    def load_config_colour_dict(self, config_key, default_colours):
        """Loads the specified colour dict from the config, resorting to the given defaults if a key is not found or is malformed."""
        colours = {}
        for cls, default_colour in default_colours.items():
            colour = self.app.config.get(config_key, cls, fallback=None)
            if not colour:
                colour = default_colour
                self.app.config.set(config_key, cls, colour)
            try:
                r, g, b, a = colour.split(',')
            except:
                print('WARNING: colour value for config key {} invalid. Resorting to default.'.format(colour))
                self.app.config.set(config_key, cls, default_colour)
                r, g, b, a = default_colour.split(',')
            colour_tup = (float(r)/255, float(g)/255, float(b)/255, float(a)/255)
            colours[cls] = colour_tup
        return colours


# # #  Main elements  # # #
class AlnKeyApp(App):
    def __init__(self, *args, **kwargs):
        # self.config - automatically added
        # self.root - automatically added; is the ScreenManager after build() completes
        # self.user_data_dir # cross-platform dir to store persistent app data
        self.config_path_keys = set() # Populated by the Save and Load base classes
        super().__init__(*args, **kwargs)
    def build(self):
        self.icon = os.path.join('app_data', 'icon_small.png')
        return Builder.load_file(os.path.join('app_data', 'aln_key_layout.kv'))
    def build_config(self, config):        
        for config_key in config_category_keys:
            config.setdefaults(config_key, {})
    def get_application_config(self):
        app_name = type(self).__name__.lower()
        if app_name.endswith('app'):
            app_name = app_name[:-3]
        return os.path.join(self.user_data_dir, '{}.ini'.format(app_name))
    def on_start(self):
        self.validate_config_paths()
        # Colours are validated when each colour dict is loaded.
    def on_stop(self):
        # may be called twice if I use app.stop(); a known bug in Kivy.
        self.config.write()
    def validate_config_paths(self):
        # In case the filesystem changes between uses of the app, config gets corrupted, etc
        default_path = os.path.expanduser("~")
        if not os.path.isdir(default_path):
            default_path = os.path.abspath('.')
        for config_key in self.config_path_keys:
            path = self.config.get(config_paths_category, config_key, fallback=None)
            if path==None or not os.path.isdir(path):
                self.config.set(config_paths_category, config_key, default_path)


class AlnKeyMain(ScreenManager):
    # From within a screen, self.manager returns this, self.manager.get_screen('name') gets other screens
    # In kivylang, accessed by `app.root`
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.app = App.get_running_app()
        self.main_screen = MainScreen(name='MainScreen')
        self.variation_screen = VariationScreen(name='VariationScreen')
        self.pdb_screen = PDBScreen(name='PDBScreen')
        self.option_screen = OptionScreen(name='OptionScreen')
        self.add_widget(self.main_screen)
        self.add_widget(self.variation_screen)
        self.add_widget(self.pdb_screen)
        self.add_widget(self.option_screen)
        #self.screens = [self.main_screen, self.variation_screen]
        # #  Common popups
        self.load_popup = LoadFilePopup()
        self.save_popup = SaveFilePopup()
        self.permission_popup = PermissionPopup()
        # #  Common data
        self.alignmant_path = ''
        self.alignment = None # A SeqList once loaded by MainScreen
        self.alignment_lengths = []
        self.alignment_consensus = ''
    def change_screen(self, new_screen_name):
        print('\n- changing to', new_screen_name, 'sizes:', Window.size, self.get_screen(new_screen_name).screen_size)
        self.current_screen.screen_size = Window.size
        if new_screen_name == self.main_screen.name:
            self.transition.direction = 'right'
        else:
            self.transition.direction = 'left'
        Window.size = self.get_screen(new_screen_name).screen_size
        self.current = new_screen_name


# # #  Screen classes  # # #
class MainScreen(BaseScreen, LoadFiles, SaveFiles):
    default_aln_status_text = 'Load an alignment file to begin.\n\nValid formats: FASTA, Clustal, & Phylip.'
    aln_status_text = StringProperty(default_aln_status_text)
    aln_missing = BooleanProperty(True)
    def __init__(self, *args, **kwargs):
        self.path_config_keys = ['load_alignment', 'save_alignment']
        self.load_alignment_key = self.path_config_keys[0]
        self.save_alignment_key = self.path_config_keys[1]
        super().__init__(*args, **kwargs)
        self.export_dropdown = AlnExportDropDown(self)
    # #  Info display
    def alignment_loaded(self):
        self.display_alignment_info()
    def display_alignment_info(self):
        aln_max_name_len = 25
        disp_fname = os.path.basename(self.manager.alignmant_path)
        if len(disp_fname) > aln_max_name_len:
            disp_fname = '...{}'.format(disp_fname[3-aln_max_name_len:])
        elif len(disp_fname) <= aln_max_name_len - 2:
            disp_fname = '  ' + disp_fname
        aln = self.manager.alignment
        num_seqs, aln_len, gap_perc = len(aln), len(aln[0]), round(aln.gaps/aln.size*100)
        aln_min, aln_max = min(self.manager.alignment_lengths), max(self.manager.alignment_lengths)
        len_str = '{}-{}'.format(aln_min, aln_max)
        self.aln_status_text = 'Loaded file:\n[b]{}[/b]\nSequences: [b]{}[/b]\nLength range: [b]{}[/b]\nAlignment length: [b]{}[/b]\nGap content: [b]{}%[/b]'.format(disp_fname, num_seqs, len_str, aln_len, gap_perc)
    # #  Loading / saving functions
    def load_alignment_button(self):
        title = 'Select protein alignment file in FASTA, Clustal, or Phylip format'
        self.get_load_filepath(title, self.load_alignment_key, self.load_alignment)
    def load_alignment(self, path, filepath):
        if self.validate_alignment(filepath):
            for screen in self.manager.screens:
                screen.alignment_loaded()
    def validate_alignment(self, filepath):
        self.manager.alignmant_path = ''
        self.manager.alignment = None
        self.manager.alignment_lengths = []
        self.manager.alignment_consensus = ''
        self.aln_missing = True
        try:
            aln = sequ.load(filepath)
        except:
            self.aln_status_text = '[b]Error:[/b] could not load chosen alignment file.\n\n{}'.format(self.default_aln_status_text)
            return False
        # Check if alignment appears valid
        if not aln.is_alignment:
            self.aln_status_text = '[b]Error:[/b] invalid alignment file due to sequences of differing lengths.\n\n{}'.format(self.default_aln_status_text)
            return False
        self.manager.alignmant_path = filepath
        self.manager.alignment = aln
        self.manager.alignment_lengths = aln.lengths
        self.manager.alignment_consensus = align.consensus(aln)
        self.aln_missing = False
        return True
    def export_alignment_button(self, format):
        title = 'Choose where to export the alignment'
        filename = os.path.basename(self.manager.alignmant_path)
        if not filename:
            filename = 'alignment'
        elif '.' in filename:
            filename = filename.rpartition('.')[0]
        export_args = ()
        # #  Format cases
        if format == 'fasta':
            filename += '.aln'
            export_fxn = sequ.save_fasta
        elif format == 'clustal':
            filename += '.clustal'
            export_fxn = sequ.save_clustal
        elif format == 'phylip interleaved':
            filename += '.phylip'
            export_fxn = sequ.save_phylip
            export_args = ('interleaved', False)
        elif format == 'phylip strict':
            filename += '.phylip'
            export_fxn = sequ.save_phylip
            export_args = ('interleaved', True)
        elif format == 'phylip sequential':
            filename += '.phylip'
            export_fxn = sequ.save_phylip
            export_args = ('sequential', False)
        elif format == 'fasta sequences':
            title = 'Choose where to export the sequences'
            filename += '.faa'
            export_fxn = sequ.save_fasta_sequences
        callback_args = (export_fxn,) + export_args
        self.get_save_filepath(title, filename, None, self.save_alignment_key, self.export_alignment, callback_args=callback_args)
    def export_alignment(self, dirpath, filename, export_fxn, *export_args):
        filepath = os.path.join(dirpath, filename)
        export_fxn(self.manager.alignment, filepath, *export_args)


class VariationScreen(BaseScreen, SaveFiles, DrawGraphics):
    focused_sequence = StringProperty('Alignment consensus') # The identity of the displayed sequence
    hover_text = StringProperty('') # Label used to display mouseover information
    def __init__(self, *args, **kwargs):
        self.path_config_keys = ['variation_save_image']
        self.save_image_key = self.path_config_keys[0]
        super().__init__(*args, **kwargs)
        self.screen_size = (1000, 800)
        self.variations = []
        self.variants = []
        self.hover_elements = {}
        default_variation_colours = {'segment_background':'255,255,255,255', 'graph_line':'17,64,67,255', 'mean_line':'150,150,150,255', 'residue_variation':'10,239,255,255', 'residue_default':'182,182,182,255', 'residue_font':'50,50,50,255', 'axis_font':'50,50,50,255'}
        self.register_colour_dict(config_variation_colours_category, default_variation_colours)
        self.sequence_selector = AlnSelectChoice(self, on_change=self.change_focused_sequence, default_choice=self.focused_sequence)
        # #  Added in aln_key_layout.kv
        #self.draw_canvas = BoxLayout
        #self.draw_canvas_view = ScrollView
        #self.show_sequence_cb = Bool
        #self.graph_type_radio = Str
        #self.num_variants_input = LabelledInput
        #self.colour_by_radio = Str
        #self.first_res_input = LabelledInput
        #self.show_range_input = LabelledInputRange
        #self.show_meanline_cb = Bool
        #self.show_numbers_cb = Bool
        #self.show_ticks_cb = Bool

        Window.bind(mouse_pos=self.canvas_mouseover) # TODO does this need to go into on_enter? or can multiple widgets be bound at once?
    
    def on_enter(self):
        if len(self.sequence_selector.container.children) == 0:
            self.sequence_selector.add(self.sequence_selector.default_choice)
            self.focused_sequence = self.sequence_selector.default_choice
            for name in self.manager.alignment.names:
                self.sequence_selector.add(name)

        if not self.variations:
            print('entered')
            self.draw_canvas.canvas.clear()
            self.draw_loading_graphic()
            
            t = threading.Thread(target=self.load_variation_then_draw)
            #t.daemon = True
            t.start()
            print('check 1')
            
            #Clock.schedule_once(self.load_variation_then_draw, 0) # Have to do this or the graphic doesn't display before the thread locks up

    def alignment_loaded(self):
        # Reset backend objects
        self.variations = []
        self.consensus = ''
        self.variants = []
        self.sequence_selector.clear()
        # Reset input objects
        self.first_res_input.reset()
        self.show_range_input.reset()
        # Reset drawing objects
        self.draw_canvas.canvas.clear()
        self.clear_hover_messages()
        
    def load_variation_then_draw(self, dt=None):
        print('loading', dt)
        self.variations = align.variation(self.manager.alignment) # long step
        self.variants = align.variants(self.manager.alignment) # brief but noticable pause
        # Can I combine both of the above functions into a single pass through the alignment?
        print('done loading')
        #self.draw_graphics()
        Clock.schedule_once(self.draw_graphics, 0) # Clock calls are performed by the main thread, and only the main thread can draw.
    
    def change_focused_sequence(self, name):
        self.focused_sequence = name
    
    # #  Mouseover methods
    def canvas_mouseover(self, window, pos):
        if not self.hover_elements or self.manager.current_screen != self.manager.variation_screen:
            return False
        if self.draw_canvas_view.collide_point(*pos):
            canv_x, canv_y = self.draw_canvas_view.to_local(*pos)
            self.hover_text = self.get_hover_message(canv_x, canv_y)
    def clear_hover_messages(self):
        self.hover_elements = {}
    def add_hover_message(self, x1, x2, y1, y2, message):
        self.hover_elements.setdefault((y1, y2), []).append((x1, x2, message))
    def process_hover_messages(self):
        """Converts the dict into a sorted list for quick lookups."""
        msgs = []
        for y1, y2 in sorted(self.hover_elements):
            row_msgs = [(y1, y2)]
            for x1, x2, msg in sorted(self.hover_elements[(y1,y2)]):
                row_msgs.append((x1, x2, msg))
            msgs.append(row_msgs)
        self.hover_elements = msgs
    def get_hover_message(self, canv_x, canv_y):
        """canv_x and canv_y are the x,y coordinates of the mouse in terms relative to the canvas."""
        #self.hover_elements = [[(y1, y2), (x1, x2, msg), ...], [(y1,y2),...]]
        for row_msgs in self.hover_elements:
            y1, y2 = row_msgs[0]
            if canv_y < y1:
                return ''
            elif canv_y > y2:
                continue
            else:
                for x1, x2, msg in row_msgs[1:]:
                    if canv_x < x1:
                        return ''
                    elif canv_x > x2:
                        continue
                    else:
                        return msg
        return ''

    # #  I/O methods
    def save_image_button(self):
        title = 'Choose where to save an image of the alignment variation'
        self.get_save_filepath(title, 'alignment_variation.png', '.png', self.save_image_key, self.save_image)
    def save_image(self, dirpath, filename):
        self.draw_canvas.export_to_png(os.path.join(dirpath, filename))

    # #  Graphical methods
    def draw_loading_graphic(self, dt=None):
        print('drawing')
        self.draw_canvas_view.scroll_x = 0
        self.draw_canvas_view.scroll_y = 1
        draw_label = self.draw_methods_factory()['draw_label']
        text = 'Calculating alignment variation...'
        label_pos = (10, self.draw_canvas_view.height-10)
        font_colour = (50/255, 50/255, 50/255, 1)
        with self.draw_canvas.canvas:
            draw_label(label_pos, text, halign='left', valign='top', font_size=26, font_colour=font_colour)
            #self.cb = Callback(self.load_variation_then_draw)
        #self.draw_canvas.canvas.ask_update()
        
    def draw_graphics(self, dt=None):
        # Get basic options
        graph_type = self.graph_type_radio
        sequence_colours = self.colour_by_radio
        first_label_num = int(self.first_res_input.text)
        # Get dynamic options
        if graph_type == 'line w variants':
            graph_type = 'line'
            num_variants = int(self.num_variants_input.text)
        else:
            num_variants = 0


        ignore_gaps = self.filter_gaps_cb # False #True # Should be set by UI element
        # TODO: Allow selection of sequence to display
        sequence_name = 'Consensus'
        #sequence_name = 'L6-M982-M982_Nme'


        if sequence_name == 'Consensus':
            raw_sequence = self.manager.alignment_consensus
        else:
            raw_sequence = self.manager.alignment.get(sequence_name).sequence
        # Format working data
        seq_ind, num_seqs, max_variants = first_label_num, len(self.manager.alignment), 0
        # sequence = [(residue, seq_ind, percent of column), ...]
        sequence, variations, variants = [], [], []
        for i in range(len(raw_sequence)):
            residue = raw_sequence[i]
            if ignore_gaps and residue == '-':
                continue
            for c, count in self.variants[i]:
                if c == residue:
                    perc = count / num_seqs * 100
                    break
            else:
                print('\n\nSomething weird happened with the variants')
            sequence.append((residue, seq_ind, perc))
            if residue != '-':
                seq_ind += 1
            variations.append(self.variations[i])
            if num_variants:
                vnts_norm = self.variants[i][0][1] # Normalized against most frequent character, not total
                if sequence_name == 'Consensus':
                    vnts = [(c, count/vnts_norm, count/num_seqs*100) for c, count in self.variants[i][1:num_variants+1]]
                else:
                    vnts = [(c, count/vnts_norm, count/num_seqs*100) for c, count in self.variants[i][:num_variants]]
                variants.append(vnts)
                max_variants = max(len(vnts), max_variants)
        # Deal with data bounds
        seq_start, seq_end = self.show_range_input.texts
        final_res = first_label_num + len(sequence) - 1
        if not seq_start:
            seq_start = first_label_num
        else:
            seq_start = int(seq_start)
            if seq_start < first_label_num:
                seq_start = first_label_num
                self.show_range_input.text1 = str(seq_start)
        if not seq_end:
            seq_end = final_res
        else:
            seq_end = int(seq_end)
            if seq_end > final_res or seq_end <= seq_start:
                seq_end = final_res
                self.show_range_input.text2 = str(seq_end)
        # Crop sequence & scores & variants
        if seq_end - seq_start + 1 < len(sequence):
            crop_start = seq_start - first_label_num
            crop_end = len(sequence) - (final_res - seq_end)
            first_label_num = seq_start
            sequence = sequence[crop_start:crop_end]
            variations = variations[crop_start:crop_end]
            variants = variants[crop_start:crop_end]

        self.draw_sequence_variation(sequence, variations, variants=variants, max_variants=max_variants, show_sequence=self.show_sequence_cb, graph_type=graph_type, sequence_colours=sequence_colours, show_mean_line=self.show_meanline_cb, show_numbers=self.show_numbers_cb, show_ticks=self.show_ticks_cb, first_label_num=first_label_num)

    def draw_sequence_variation(self, sequence, variations, variants=[], max_variants=0, show_sequence=True, graph_type='line', sequence_colours='properties', show_mean_line=True, show_numbers=True, show_ticks=True, first_label_num=1):
        """Auto-adapts to the current screen size as resized by the user."""
        # #  Parameters
        x_padding, y_padding = 3, 3 # Space between graphics and edge
        segment_spacing = 10 # Vertical space between each segment
        segment_x_padding, segment_y_padding = 2, 2 # Space between segment edge and res/graph
        graph_buffer = 5 # Space between residues and graph
        dash_len_space = (10, 10) # Dash line length, space between dashes
        inds_tick = 5
        inds_buffer = 1 # Space between residues and indices
        inds_char_width, inds_height = 0, 0 # Calculated if needed
        colours = self.get_colour_dict(config_variation_colours_category)
        res_font_clr = colours['residue_font']
        ax_font_clr = colours['axis_font']
        default_res_bkgrnd = colours['residue_default']

        # #  Compute drawing parameters
        aln_len = len(sequence)
        max_var = max(variations)
        canv_w, canv_h = self.draw_canvas_view.size
        if show_sequence:
            res_size = 20, 20
            res_spacing = 1 # Space between residues
            minor_tick_interval = 5
            major_tick_interval = 10
            res_per_segment = int((canv_w-2*x_padding-2*segment_x_padding+res_spacing)//(res_size[0]+res_spacing))
            num_segs = aln_len // res_per_segment
            if aln_len % res_per_segment > 0:
                num_segs += 1
        else:
            avail_seg_w = canv_w - 2*x_padding - 2*segment_x_padding
            res_w = max(avail_seg_w / aln_len, 2)
            res_size = res_w, 20
            res_spacing = 0
            major_tick_interval = pick_nice_interval(aln_len)
            minor_tick_interval = major_tick_interval / 2
            res_per_segment = aln_len
            num_segs = 1
        if variants and max_variants > 0:
            graph_height = max_variants * res_size[1]
        else:
            graph_height = 60
        if show_mean_line:
            mean_val = sum(variations) / len(variations)
            mean_height = mean_val/max_var * graph_height
        if show_numbers:
            test_lab = CoreLabel(text='0', bold=True)
            test_lab.refresh()
            inds_char_width, inds_height = test_lab.texture.size
            del test_lab
        if not show_ticks:
            inds_tick, inds_buffer = 0, 0
        segment_width = 2*segment_x_padding + res_per_segment*res_size[0] + (res_per_segment-1)*res_spacing
        segment_height = 2*segment_y_padding + inds_height + inds_tick + inds_buffer + res_size[1] + graph_buffer + graph_height
        new_canv_w = segment_width + 2*x_padding
        new_canv_h = num_segs*segment_height + 2*y_padding
        if num_segs > 1:
            new_canv_h += (num_segs-1) * segment_spacing

        # #  Do drawing
        self.clear_hover_messages()
        self.draw_canvas.canvas.clear()
        self.draw_canvas.width = new_canv_w
        self.draw_canvas.height = new_canv_h
        draw_methods = self.draw_methods_factory()
        draw_line, draw_curve, draw_rect, draw_label = draw_methods['draw_line'], draw_methods['draw_curve'], draw_methods['draw_rect'], draw_methods['draw_label']
        with self.draw_canvas.canvas:
            seg_x = x_padding
            seg_left = seg_x + segment_x_padding
            line_coords, prev_graph_h = [], None
            for seg_num in range(num_segs):
                # Segment values
                seg_start = seg_num * res_per_segment
                seg_end = min(seg_start + res_per_segment, aln_len)
                num_res = seg_end - seg_start
                seg_y = new_canv_h-y_padding-segment_height-(segment_height+segment_spacing)*seg_num
                seg_w = 2*segment_x_padding + num_res*res_size[0] + (num_res-1)*res_spacing
                seg_right = seg_x + seg_w-segment_x_padding
                # Other element values
                res_y = seg_y + segment_y_padding + inds_height + inds_tick + inds_buffer
                graph_y = res_y + res_size[1] + graph_buffer
                # Draw background
                draw_rect((seg_x, seg_y), (seg_w, segment_height), colours['segment_background'])
                prev_seg_pnt = None
                ind_coords, seg_line_coords = [], []
                for seg_ind in range(seg_end - seg_start):
                    ind = seg_start + seg_ind
                    residue = sequence[ind][0]
                    var_prop = variations[ind]/max_var
                    # Draw residue
                    res_x = seg_left + (res_size[0]+res_spacing)*seg_ind
                    res_x_mid = res_x + res_size[0]/2
                    if sequence_colours == 'properties':
                        res_clr = self.get_residue_colour(residue)
                    elif sequence_colours == 'conservation':
                        res_clr = colours['residue_variation'][:3] + (1-var_prop,)
                    else:
                        res_clr = None
                    res_name = residue if residue != '-' else 'Gap prior to '
                    hover_message = '{}{} ({}%)'.format(res_name, sequence[ind][1], round(sequence[ind][2]))
                    if show_sequence:
                        rect = draw_label((res_x,res_y), residue, size=res_size, font_colour=res_font_clr, box_colour=res_clr)
                        self.add_hover_message(res_x, res_x+res_size[0], res_y, res_y+res_size[1], hover_message)
                    elif res_clr:
                        rect = draw_rect((res_x,res_y), res_size, res_clr)
                        self.add_hover_message(res_x, res_x+res_size[0], res_y, res_y+res_size[1], hover_message)
                    # Deal with graphs
                    graph_h = var_prop * graph_height
                    if graph_type == 'bar':
                        draw_rect((res_x,graph_y), (res_size[0],graph_h), colours['graph_line'])
                    elif graph_type == 'line':
                        if seg_ind == 0: # First residue of the segment
                            if prev_graph_h == None: # First segment
                                prev_graph_y = graph_y + graph_h
                            else:
                                prev_graph_y = graph_y + prev_graph_h
                            prev_seg_pnt = res_x_mid-res_spacing-res_size[0]/2, prev_graph_y
                        seg_line_coords.append((res_x_mid, graph_y+graph_h))
                        prev_graph_h = graph_h
                        # Graph variant residues
                        if variants:
                            vnt = variants[ind]
                            vnt_ind = 0
                            for c_res, c_prop, c_perc in vnt:
                                c_res_y = graph_y + vnt_ind*res_size[1]
                                if sequence_colours == 'properties':
                                    c_res_clr = self.get_residue_colour(c_res) or default_res_bkgrnd
                                else:
                                    c_res_clr = default_res_bkgrnd
                                res_op = max(c_prop, 0.05)
                                c_res_clr = c_res_clr[:3] + (res_op,)
                                if show_sequence:
                                    draw_label((res_x,c_res_y), c_res, size=res_size, font_colour=res_font_clr, box_colour=c_res_clr)
                                else:
                                    draw_rect((res_x,c_res_y), res_size, c_res_clr)
                                hover_message = '{} ({}%)'.format(c_res if c_res != '-' else 'Gap', round(c_perc))
                                self.add_hover_message(res_x, res_x+res_size[0], c_res_y, c_res_y+res_size[1], hover_message)
                                vnt_ind += 1
                    # Record values for indices
                    res_num = first_label_num + ind
                    if res_num % major_tick_interval == 0:
                        ind_coords.append((res_x_mid, res_y-inds_buffer, res_num))
                    elif res_num % minor_tick_interval == 0:
                        ind_coords.append((res_x_mid, res_y-inds_buffer, None))
                # Deal with segment line graph start/ends
                if graph_type == 'line':
                    if ind < aln_len - 1:
                        next_pnt_var = variations[ind+1]
                    else: # For the very last point
                        next_pnt_var = variations[ind]
                    next_seg_h = next_pnt_var/max_var * graph_height
                    next_seg_pnt = res_x+res_size[0]+res_spacing, next_seg_h+graph_y
                    line_coords.append((prev_seg_pnt, next_seg_pnt, seg_line_coords))
                # Draw indices
                if show_ticks or show_numbers:
                    for tick_x, tick_y, num in ind_coords:
                        if show_ticks and (show_numbers or num != None):
                            tick_points = [(tick_x,tick_y), (tick_x,tick_y-inds_tick)]
                            draw_line(tick_points, colour=ax_font_clr)
                        if show_numbers and num != None:
                            num_str = str(num)
                            half_w = inds_char_width*len(num_str)/2
                            inds_x = max(seg_left+half_w, tick_x)
                            inds_x = min(seg_right-half_w, inds_x)
                            inds_y = tick_y - inds_tick
                            draw_label((inds_x,inds_y), num_str, halign='middle', valign='top', font_colour=ax_font_clr)
                # Draw mean line
                if show_mean_line:
                    meanline_y = mean_height + graph_y
                    meanline_points = [(seg_left+1, meanline_y), (seg_right, meanline_y)]
                    draw_line(meanline_points, width=1.0, colour=colours['mean_line'], dashes=dash_len_space)
                # Draw axes
                ax_points = [(seg_left+1, graph_y+graph_height), (seg_left+1, graph_y), (seg_right, graph_y)]
                draw_line(ax_points, width=1.1, colour=res_font_clr)
            # Draw line graph over whole segment
            if graph_type == 'line':
                for prev_seg_pnt, next_seg_pnt, seg_line_coords in line_coords:
                    draw_curve(seg_line_coords, width=2.0, colour=colours['graph_line'], continue_from=prev_seg_pnt, continue_to=next_seg_pnt)
        self.process_hover_messages()


class PDBScreen(BaseScreen):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.screen_size = (1000, 800)


class OptionScreen(BaseScreen):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


if __name__ == '__main__':
    if hasattr(sys, '_MEIPASS'):
        resource_add_path(os.path.join(sys._MEIPASS))
    AlnKeyApp().run()
