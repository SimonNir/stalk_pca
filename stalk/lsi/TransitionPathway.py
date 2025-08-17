#!/usr/bin/env python3
'''TransitionStateSearch class for finding transition pathways.'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from numpy import sort, dot, array
from stalk.lsi.PathwayImage import PathwayImage
from stalk.params.ParameterSet import ParameterSet
from stalk.util.util import directorize

class TransitionPathway():
    _images: list[PathwayImage] = []  # list of active PathwayImage objects
    _all_images: list[PathwayImage] = [] # Full set of NEB images
    _active_indices: list[int] = [] # Indices of images being optimized (do not use negatives)
    _path = ''  # base path

    def __init__(self, path='', all_images: list = None, active_indices: list[int] = None):
        self._path = path

        # Wrap only ParameterSet objects; leave PathwayImage as-is
        self._all_images = []
        if all_images:
            for img in all_images:
                if isinstance(img, ParameterSet):
                    self._all_images.append(PathwayImage(img))
                elif isinstance(img, PathwayImage):
                    self._all_images.append(img)
                else:
                    raise TypeError(f"All images must be ParameterSet or PathwayImage, got {type(img)}")
        
        self._images = []

        # Determine highest-energy image
        # Determine highest-energy image, allowing None for endpoints only
        energies = []
        for i, img in enumerate(self._all_images):
            energy = img.surrogate_energy
            if energy is None and i not in (0, len(self._all_images) - 1):
                raise ValueError(f"Image at index {i} has surrogate_energy=None (only endpoints may be None)")
            energies.append(energy if energy is not None else float('-inf')) # ignore if endpoints are None
        self._max_energy_index = max(range(len(self._all_images)), key=lambda i: energies[i])

        # Determine active indices
        if active_indices is None:
            print("No active indices specified; using eqm1, eqm2, and saddle only by default")
            self._active_indices = [0, self._max_energy_index, len(self._all_images)-1]
        else:
            self._active_indices = active_indices
            required = {0, self._max_energy_index, len(self._all_images)-1}
            if not required.issubset(self._active_indices):
                print(f"Warning: Active indices do not include first, last, and highest-energy images: {required}")

        # Add images to self._images
        for i in self._active_indices:
            image = self._all_images[i]
            if i == 0 or i == len(self._all_images)-1:
                img_type = 'eqm'
            elif i == self._max_energy_index:
                img_type = 'saddle'
            else:
                img_type = 'intermediate'
            image.image_type = img_type
            self._images.append(image)

        # calculate tangents 
        self.calculate_tangents()

    @property
    def path(self):
        return self._path

    @path.setter
    def path(self, path):
        if isinstance(path, str):
            self._path = directorize(path)
        else:
            raise TypeError("path must be a string")

    @property
    def images(self):
        return self._images
    
    @property
    def all_images(self):
        return self._all_images

    @property
    def intermediate_images(self):
        return self._images[1:-1]

    @property
    def pointA(self):
        return self._images[0] if len(self) >= 1 else None

    @property
    def pointB(self):
        return self._images[-1] if len(self) >= 2 else None

    @property
    def saddle(self):
        return self._all_images[self._max_energy_index]

    @property
    def highest_energy_active_image(self):
        return max(self.images, key=lambda i: i.structure.surrogate_energy)

    @property  # DEPRECATED
    def difference(self):
        if self.pointB is not None:
            return self.pointB.structure.params - self.pointA.structure.params
        return None

    @property
    def pathway_init(self):
        params = []
        params_err = []
        for image in self.images:
            params.append(image.structure_init.params)
            params_err.append(image.structure_init.params_err)
        return array(params), array(params_err)

    @property
    def pathway_final(self):
        params = []
        params_err = []
        for image in self.images:
            params.append(image.structure_final.params)
            params_err.append(image.structure_final.params_err)
        return array(params), array(params_err)

    def calculate_tangents(self):
        """Compute tangents for all active images using robust local energy rules, 
        originally developed by Henkelman and Jonsson, JCP, 2000."""
        n = len(self._all_images)
        for idx, img in enumerate(self._images):
            # Find full index in all_images
            full_idx = self._active_indices[idx]
            # Skip endpoints
            if full_idx == 0:
                # Forward difference for first image
                tangent = self._all_images[full_idx + 1].structure.params - img.structure.params
            elif full_idx == n - 1:
                # Backward difference for last image
                tangent = img.structure.params - self._all_images[full_idx - 1].structure.params
            else:
                # Intermediate image
                prev_img = self._all_images[full_idx - 1]
                next_img = self._all_images[full_idx + 1]
                curr_img = img

                E_prev = prev_img.surrogate_energy
                E_curr = curr_img.surrogate_energy
                E_next = next_img.surrogate_energy

                P_prev = prev_img.structure.params
                P_curr = curr_img.structure.params
                P_next = next_img.structure.params

                # First intermediate: if no energy for eqm1, assume E_curr > E_prev
                if full_idx == 1:
                    if E_prev is None:
                        E_prev = E_curr
                    tangent = P_next - P_curr
                # Last intermediate: if no energy for eqm2, assume E_curr > E_next
                elif full_idx == n - 2:
                    if E_next is None:
                        E_next = E_curr
                    tangent = P_curr - P_prev
                elif E_next > E_curr and E_prev < E_curr:
                    tangent = P_next - P_curr
                elif E_prev > E_curr and E_next < E_curr:
                    tangent = P_curr - P_prev
                elif E_next > E_curr and E_prev > E_curr:
                    # Local minimum: weighted average, weights = energy rises
                    w_fwd = E_next - E_curr
                    w_bwd = E_prev - E_curr
                    tangent = (w_fwd * (P_next - P_curr) + w_bwd * (P_curr - P_prev)) / (w_fwd + w_bwd)
                elif E_next < E_curr and E_prev < E_curr:
                    # Local maximum: weighted average, weights = energy drops
                    w_fwd = E_curr - E_next
                    w_bwd = E_curr - E_prev
                    tangent = (w_fwd * (P_next - P_curr) + w_bwd * (P_curr - P_prev)) / (w_fwd + w_bwd)
                else:
                    # Fallback: central difference
                    tangent = 0.5 * (P_next - P_prev)
            img.tangent = tangent

    def calculate_hessians(
        self,
        **hessian_args,
    ):
        for i, image in enumerate(self.images):
            image.calculate_hessian(
                tangent=image.tangent,
                path=('{}image_{}/hessian/').format(self.path, self._active_indices[i]),
                **hessian_args
            )
        # end for
    # end def

    def generate_surrogates(self, **surrogate_args):
        for image in self.images:
            image.generate_surrogate(**surrogate_args)

    def optimize_surrogates(self, **optimize_args):
        for image in self.images:
            image.optimize_surrogate(**optimize_args)

    def run_linesearches(self, **lsi_args):
        for image in self.images:
            image.run_linesearch(**lsi_args)

    # Deprecated
    def _calculate_rc(self, image: ParameterSet):
        rc = (dot(self.difference, image.params - self.pointA.structure.params) /
              dot(self.difference, self.difference))
        return rc

    def __len__(self):
        return len(self.images)
