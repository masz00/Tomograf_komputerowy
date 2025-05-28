import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import pydicom
from pydicom.dataset import Dataset, FileDataset
from pydicom.uid import ExplicitVRLittleEndian, CTImageStorage
from skimage.transform import radon, iradon
from skimage.util import img_as_ubyte
from skimage.exposure import rescale_intensity
from skimage.io import imread




def bresenham_line(x0, y0, x1, y1):
    points_x = []
    points_y = []
    
    dx = abs(x1 - x0)
    dy = abs(y1 - y0)
    sx = 1 if x0 < x1 else -1
    sy = 1 if y0 < y1 else -1
    err = dx - dy
    
    while True:
        points_x.append(x0)
        points_y.append(y0)
        
        if x0 == x1 and y0 == y1:
            break
            
        e2 = 2 * err
        if e2 > -dy:
            err -= dy
            x0 += sx
        if e2 < dx:
            err += dx
            y0 += sy
            
    return np.array(points_y), np.array(points_x)

import numpy as np
from scipy.ndimage import rotate
from scipy.fftpack import fft, ifft, fftfreq

def generate_parallel_rays(img_shape, angle_deg, detector_count, detector_span_deg):
    angle_rad = np.deg2rad(angle_deg)
    detector_span_rad = np.deg2rad(detector_span_deg)
    
    center = np.array(img_shape) / 2
    max_length = max(img_shape)

    detector_positions = np.linspace(-0.5, 0.5, detector_count) * detector_span_rad
    ray_vectors = np.array([
        [np.cos(angle_rad + np.pi / 2 + offset), np.sin(angle_rad + np.pi / 2 + offset)]
        for offset in detector_positions
    ])

    starts = []
    ends = []
    for vec in ray_vectors:
        midpoint = center + max_length * 0.5 * np.array([np.cos(angle_rad), np.sin(angle_rad)])
        start = midpoint + vec * max_length
        end = midpoint - vec * max_length
        starts.append(start.astype(int))
        ends.append(end.astype(int))
    return starts, ends

def radon_transform_parallel(img, angles, detector_count=180, detector_span=180):
    img = img.astype(np.float32)
    sinogram = np.zeros((detector_count, len(angles)), dtype=np.float32)
    img_center = np.array(img.shape) / 2
    diag = int(np.ceil(np.hypot(*img.shape)))
    max_offset = diag // 2

    for i, angle in enumerate(angles):
        rad = np.deg2rad(angle)
        sin_a, cos_a = np.sin(rad), np.cos(rad)

        for det_idx in range(detector_count):
            offset = (det_idx - detector_count / 2) * (max_offset * 2 / detector_count)
            x0 = img_center[1] + offset * cos_a - diag * 0.5 * sin_a
            y0 = img_center[0] - offset * sin_a - diag * 0.5 * cos_a
            x1 = img_center[1] + offset * cos_a + diag * 0.5 * sin_a
            y1 = img_center[0] - offset * sin_a + diag * 0.5 * cos_a

            rr, cc = bresenham_line(int(y0), int(x0), int(y1), int(x1))
            valid = (rr >= 0) & (rr < img.shape[0]) & (cc >= 0) & (cc < img.shape[1])
            sinogram[det_idx, i] = np.sum(img[rr[valid], cc[valid]])

    return sinogram

def ramp_filter(sinogram):
    projections = sinogram.T
    num_detectors = projections.shape[1]
    freqs = fftfreq(num_detectors).reshape(1, -1)
    filter_kernel = 2 * np.abs(freqs)
    projections_fft = fft(projections, axis=1)
    projections_filtered = np.real(ifft(projections_fft * filter_kernel, axis=1))
    return projections_filtered.T

def inverse_radon_transform_parallel(sinogram, angles, output_size):
    filtered_sinogram = ramp_filter(sinogram)
    reconstruction = np.zeros((output_size, output_size), dtype=np.float32)
    img_center = np.array([output_size / 2, output_size / 2])
    diag = int(np.ceil(np.hypot(output_size, output_size)))
    max_offset = diag // 2
    detector_count = sinogram.shape[0]

    norm_map = np.zeros_like(reconstruction)

    for i, angle in enumerate(angles):
        rad = np.deg2rad(angle)
        sin_a, cos_a = np.sin(rad), np.cos(rad)

        for det_idx in range(detector_count):
            offset = (det_idx - detector_count / 2) * (max_offset * 2 / detector_count)
            x0 = img_center[1] + offset * cos_a - diag * 0.5 * sin_a
            y0 = img_center[0] - offset * sin_a - diag * 0.5 * cos_a
            x1 = img_center[1] + offset * cos_a + diag * 0.5 * sin_a
            y1 = img_center[0] - offset * sin_a + diag * 0.5 * cos_a

            rr, cc = bresenham_line(int(y0), int(x0), int(y1), int(x1))
            valid = (rr >= 0) & (rr < output_size) & (cc >= 0) & (cc < output_size)
            rr = rr[valid]
            cc = cc[valid]
            reconstruction[rr, cc] += filtered_sinogram[det_idx, i]
            norm_map[rr, cc] += 1

    norm_map[norm_map == 0] = 1
    reconstruction /= norm_map
    reconstruction -= reconstruction.min()
    reconstruction /= reconstruction.max()
    return reconstruction

def convert_image_to_ubyte(img):
    return img_as_ubyte(rescale_intensity(img, out_range=(0.0, 1.0)))

def save_as_dicom(file_name, img, patient_data):
    img_converted = convert_image_to_ubyte(img)

    meta = Dataset()
    meta.MediaStorageSOPClassUID = CTImageStorage
    meta.MediaStorageSOPInstanceUID = pydicom.uid.generate_uid()
    meta.TransferSyntaxUID = ExplicitVRLittleEndian  

    ds = FileDataset(None, {}, preamble=b"\0" * 128)
    ds.file_meta = meta

    ds.is_little_endian = True
    ds.is_implicit_VR = False

    ds.SOPClassUID = CTImageStorage
    ds.SOPInstanceUID = meta.MediaStorageSOPInstanceUID

    ds.PatientName = patient_data["PatientName"]
    ds.StudyDate = patient_data["StudyDate"]
    ds.ImageComments = patient_data["ImageComments"]
    
    ds.Modality = "CT"
    ds.SeriesInstanceUID = pydicom.uid.generate_uid()
    ds.StudyInstanceUID = pydicom.uid.generate_uid()
    ds.FrameOfReferenceUID = pydicom.uid.generate_uid()

    ds.BitsStored = 8
    ds.BitsAllocated = 8
    ds.SamplesPerPixel = 1
    ds.HighBit = 7

    ds.ImagesInAcquisition = 1
    ds.InstanceNumber = 1

    ds.Rows, ds.Columns = img_converted.shape
    ds.ImageType = r"ORIGINAL\PRIMARY\AXIAL"

    ds.PhotometricInterpretation = "MONOCHROME2"
    ds.PixelRepresentation = 0

    pydicom.dataset.validate_file_meta(ds.file_meta, enforce_standard=True)
    ds.PixelData = img_converted.tobytes()

    ds.save_as(file_name, write_like_original=False)



st.title("Symulator Tomografii Komputerowej")
tab1, tab2 = st.tabs(["Symulacja TK", "Odczyt DICOM"])

with tab1:
    uploaded_file = st.file_uploader("Wybierz obraz (skala szarości)", type=["png", "jpg", "jpeg"], key="image_uploader")

    if uploaded_file:
        img = imread(uploaded_file, as_gray=True)

        angle_step = st.slider("Krok kąta (∆α)", min_value=0.5, max_value=4.0, step = 0.5, value=2.0)
        detector_count = st.slider("Liczba detektorów", min_value=90, max_value=720, step=90, value=180)
        fan_angle = st.slider("Rozpiętość wachlarza (stopnie)", min_value=45, max_value=270, step=45, value=180)
        
        theta = np.arange(0, 180, angle_step)  # angle_step z suwaka
        sinogram = radon_transform_parallel(img, theta, detector_count=detector_count, detector_span=fan_angle)
        reconstructed_img = inverse_radon_transform_parallel(sinogram, theta, output_size=img.shape[0])

        st.subheader("Obraz wejściowy")
        st.image(img, clamp=True, caption="Obraz wejściowy")

        st.subheader("Sinogram")
        fig, ax = plt.subplots()
        ax.imshow(sinogram, cmap="gray", aspect="auto")
        st.pyplot(fig)

        st.subheader("Odtworzony obraz")
        fig, ax = plt.subplots()
        ax.imshow(reconstructed_img, cmap="gray")
        st.pyplot(fig)
        
        st.subheader("Dane pacjenta")
        patient_name = st.text_input("Imię i nazwisko pacjenta", "Jan Kowalski", key="patient_name")
        study_date = st.date_input("Data badania", key="study_date")
        image_comments = st.text_area("Komentarze", "Symulacja TK", key="image_comments")

        if st.button("Zapisz jako DICOM"):
            patient_data = {
                "PatientName": patient_name,
                "StudyDate": study_date.strftime("%Y%m%d"),
                "ImageComments": image_comments
            }

            save_as_dicom(f"{patient_name.replace(' ', '_')}_{study_date}.dcm", reconstructed_img, patient_data)
            st.success(f"Plik DICOM został zapisany jako `{patient_name.replace(' ', '_')}_{study_date}.dcm`")

with tab2:
    st.subheader("Odczyt plików DICOM")
    uploaded_dicom = st.file_uploader("Wybierz plik DICOM do odczytu", type=["dcm"], key="dicom_uploader")
    
    if uploaded_dicom:
        try:
            dicom_data = pydicom.dcmread(uploaded_dicom, force=True)
            
            st.subheader("Informacje z pliku DICOM")
            col1, col2 = st.columns(2)

            def format_date(date_string):
                if date_string and len(date_string) == 8:
                    return f"{date_string[6:8]}/{date_string[4:6]}/{date_string[0:4]}"
                return "Brak danych"
            with col1:
                st.markdown("**Imię i nazwisko:**")
                st.markdown("**Data badania:**")
                st.markdown("**Komentarze:**")
            with col2:
                st.markdown(f"**{dicom_data.get('PatientName', 'Brak danych')}**")
                st.markdown(f"**{format_date(dicom_data.get('StudyDate', ''))}**")
                st.markdown(f"**{dicom_data.get('ImageComments', 'Brak danych')}**")
            if hasattr(dicom_data, 'pixel_array'):
                st.subheader("Obraz medyczny")
                st.image(dicom_data.pixel_array, caption=f"{dicom_data.get('Modality', '')} {dicom_data.get('BodyPartExamined', '')}", use_container_width=True)
            
        except Exception as e:
            st.error(f"Błąd odczytu pliku: {e}")
